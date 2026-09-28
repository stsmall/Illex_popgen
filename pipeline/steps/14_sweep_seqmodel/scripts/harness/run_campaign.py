#!/usr/bin/env python
"""Campaign runner: fan out generate.py sims + aggregate per-class training ms.

generate.py APPENDS to per-class files, so concurrent sims into one dir would
interleave. This runs each sim in its OWN dir (parallel, thread-capped), then
aggregates the per-sim per-position ms into the combined per-class training files
diploSHIC makeTrainingSets consumes:

    neutral  -> neutral.msOut.gz                 (all 11 windows/sim, all neutral)
    sweep C  -> C_0.msOut.gz .. C_10.msOut.gz    (window i = sweep in subwin i;
                                                   _5 = focal sweep, others = linked)

Usage:
    run_campaign.py --klass hard --n 60 --concurrency 16 --outdir DIR [--s S] [--F0 F]
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse       # noqa: E402
import glob           # noqa: E402
import gzip           # noqa: E402
import shutil         # noqa: E402
import subprocess     # noqa: E402
import sys            # noqa: E402
import time           # noqa: E402
from concurrent.futures import ThreadPoolExecutor, as_completed  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_GEN = os.path.join(_HERE, "generate.py")
_PY = sys.executable       # the slim_sim python running this


def run_one(klass, seed, s, F0, simdir, min_maf=None, keep_slim_trees=True):
    cmd = [_PY, _GEN, "--klass", klass, "--seed", str(seed), "--outdir", simdir]
    if s is not None:
        cmd += ["--s", str(s)]
    if F0 is not None:
        cmd += ["--F0", str(F0)]
    if min_maf is not None:
        cmd += ["--min-maf", str(min_maf)]
    cmd += ["--keep-slim-trees"] if keep_slim_trees else ["--no-keep-slim-trees"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return seed, r.returncode, (r.stderr[-400:] if r.returncode else "")


_COMBINE_SH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "combine_ms.sh")


def combine(files, out, n_hap, seed=0):
    """Combine gzipped ms files into one: a single header + every // block.

    Delegates to combine_ms.sh (zcat/gzip, C-speed). Python's gzip compresses
    at ~5 MB/s, which made this crawl ~1 h/position on the 240 MB/position pilot
    files and ran unmonitored for hours; the shell path is ~10-20x faster.
    Returns the total block count.
    """
    r = subprocess.run(["bash", _COMBINE_SH, out, str(n_hap), *files],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"combine_ms.sh failed ({r.returncode}): {r.stderr[-800:]}")
    return int(r.stdout.strip() or 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--klass", required=True,
                    choices=["neutral", "hard", "soft", "partialHard", "partialSoft"])
    ap.add_argument("--n", type=int, required=True, help="number of sims (each -> 11 windows).")
    ap.add_argument("--concurrency", type=int, default=8)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--s", type=float, default=None)
    ap.add_argument("--F0", type=float, default=None)
    ap.add_argument("--seed0", type=int, default=1000, help="first seed (seeds = seed0..seed0+n-1).")
    ap.add_argument("--n-dip", type=int, default=350, dest="n_dip")
    ap.add_argument("--min-maf", type=float, default=None, dest="min_maf",
                    help="passed to generate.py (0 -> full-SFS training set; default: config min_maf 0.01).")
    ap.add_argument("--keep-slim-trees", action=argparse.BooleanOptionalAction, default=True,
                    help="passed to generate.py; also rescues each sim's *_recon.trees into "
                         "<outdir>/trees/<klass>/ before per-sim cleanup. Default ON.")
    ap.add_argument("--keep-sims", action="store_true",
                    help="keep the per-sim dirs (default: delete after a clean aggregate to save disk).")
    a = ap.parse_args()

    simroot = os.path.join(a.outdir, "sims", a.klass)
    os.makedirs(simroot, exist_ok=True)
    seeds = list(range(a.seed0, a.seed0 + a.n))
    t0 = time.time()
    ok, fails = 0, []
    with ThreadPoolExecutor(max_workers=a.concurrency) as ex:
        futs = {ex.submit(run_one, a.klass, sd, a.s, a.F0, os.path.join(simroot, f"s{sd}"),
                          a.min_maf, a.keep_slim_trees): sd
                for sd in seeds}
        for fu in as_completed(futs):
            sd, rc, err = fu.result()
            if rc == 0:
                ok += 1
            else:
                fails.append((sd, err))
            done = ok + len(fails)
            if done % 10 == 0 or done == a.n:
                print(f"[campaign {a.klass}] {done}/{a.n} ({ok} ok, {len(fails)} fail), "
                      f"{(time.time() - t0) / 60:.1f} min", flush=True)
    if fails:
        print(f"[campaign {a.klass}] {len(fails)} FAILURES (showing 5):", flush=True)
        for sd, err in fails[:5]:
            print(f"  seed {sd}: {err}", flush=True)

    # ---- aggregate per-sim -> combined per-class training files -------------
    n_hap = 2 * a.n_dip
    if a.klass == "neutral":
        targets = {"neutral.msOut.gz": os.path.join(simroot, "s*", "neutral.msOut.gz")}
    else:
        targets = {f"{a.klass}_{i}.msOut.gz": os.path.join(simroot, "s*", f"{a.klass}_{i}.msOut.gz")
                   for i in range(11)}
    total_blocks = 0
    for name, pat in targets.items():
        files = sorted(glob.glob(pat))
        nb = combine(files, os.path.join(a.outdir, name), n_hap)
        total_blocks += nb
        print(f"[campaign] {name}: {nb} blocks from {len(files)} sims", flush=True)

    if a.keep_slim_trees and not fails:        # rescue reconstructed trees before per-sim cleanup
        treedir = os.path.join(a.outdir, "trees", a.klass); os.makedirs(treedir, exist_ok=True)
        moved = 0
        for tf in glob.glob(os.path.join(simroot, "s*", "*_recon.trees")):
            shutil.move(tf, os.path.join(treedir, os.path.basename(tf))); moved += 1
        print(f"[campaign] preserved {moved} recon.trees -> {treedir}", flush=True)
    if not a.keep_sims and not fails:
        shutil.rmtree(simroot, ignore_errors=True)
        print(f"[campaign] removed per-sim dirs ({simroot})", flush=True)
    print(f"CAMPAIGN_DONE klass={a.klass} sims_ok={ok}/{a.n} total_blocks={total_blocks} "
          f"wall={(time.time() - t0) / 60:.1f}min", flush=True)


if __name__ == "__main__":
    main()
