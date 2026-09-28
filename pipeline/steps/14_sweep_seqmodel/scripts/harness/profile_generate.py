#!/usr/bin/env python
"""Decompose the generate.py per-sim cost into stages (uncontended timing).

Answers two Phase-B budget questions: (1) how much of the ~35 min/sim measured at
5-way parallel was contention (run this alone), and (2) which stage dominates
(SLiM forward vs pyslim.recapitate+msprime overlay vs genotype matrix vs the 11x
ascertain) -- i.e. what to optimize. Reuses generate.py's helpers so it profiles
the real pipeline.

    OMP_NUM_THREADS=1 python profile_generate.py --klass hard --seed 99 --outdir DIR
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse   # noqa: E402
import sys        # noqa: E402
import time       # noqa: E402

import numpy as np  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_SLIM = os.path.join(_HERE, "slim")
for _p in (_HERE, _SLIM):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import ms_export          # noqa: E402
import ascertain          # noqa: E402
import run_slim_sweep as rss                       # noqa: E402
from neutral_msprime import simulate_neutral       # noqa: E402
from generate import _build_slim_args, KLASS_MAP   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--klass", default="hard", choices=list(KLASS_MAP))
    ap.add_argument("--seed", type=int, default=99)
    ap.add_argument("--s", type=float, default=0.01)
    ap.add_argument("--F0", type=float, default=None)
    ap.add_argument("--L", type=int, default=2_200_000)
    ap.add_argument("--n-dip", type=int, default=350, dest="n_dip")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    sweep_type, mode = KLASS_MAP[a.klass]
    T = {}

    t = time.time()
    if a.klass == "neutral":
        ts = simulate_neutral(a.seed, L=a.L, n_dip=a.n_dip)
        T["1_sim(msprime neutral)"] = time.time() - t
    else:
        args = _build_slim_args(a.seed, sweep_type, mode, a.s, a.F0, a.L, a.n_dip)
        raw = os.path.join(a.outdir, "raw.trees")
        save = os.path.join(a.outdir, "save.trees")
        est = os.path.join(a.outdir, "est.trees")
        info = rss.run_slim(args, raw, save, est_trees=est)
        T["1_slim_forward"] = time.time() - t
        print("  slim: nrestart=%s status=%s" % (info.get("nrestart"),
              info.get("status") or info.get("present_freq")), flush=True)
        t = time.time()
        ts = rss.reconstruct(args, raw)
        T["2_reconstruct(subsample+recap+overlay)"] = time.time() - t

    t = time.time()
    G, pos = ms_export.genotype_and_positions(ts)
    T["3_genotype_matrix"] = time.time() - t
    print("  n_sites=%d n_hap=%d" % (ts.num_sites, ts.num_samples), flush=True)

    t = time.time()
    xpos = a.L // 2
    keep = np.abs(pos - xpos) > 1.0
    Gf, posf = G[keep, :], pos[keep]
    model = ascertain.load_model()
    rng = np.random.default_rng(a.seed)
    for sw, off in ms_export.window_offsets(sweep_pos=xpos):
        geno, pr = ms_export.window_from(Gf, posf, off, 1_100_000)
        ascertain.ascertain(geno, pr, rng, model)
    T["4_window+ascertain_11x"] = time.time() - t

    tot = sum(T.values())
    print("=== STAGE TIMING (uncontended) klass=%s ===" % a.klass, flush=True)
    for k, v in T.items():
        print("  %-42s %8.1f s  (%4.1f%%)" % (k, v, 100 * v / tot), flush=True)
    print("  %-42s %8.1f s  (%.1f min)" % ("TOTAL", tot, tot / 60), flush=True)


if __name__ == "__main__":
    main()
