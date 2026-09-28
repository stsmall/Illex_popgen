#!/usr/bin/env python3
"""run_pilot_campaign.py -- PLAN_02 Phase A+B pilot: generate contiguous calibration segments
(SLiM hard/soft sweeps at known XPOS + msprime neutral) under the variable ReLERNN map on chr43,
then scan each through the IDENTICAL empirical fvecVcf path (calib_scan_segment.py, illexModel_vcf)
to per-window class probs. Writes a manifest the fit step consumes.

Each worker runs gen -> scan sequentially for one segment (~1 core). Concurrency-capped for the
shared box. Deterministic: sel/loci drawn from a seeded RNG; every sim seeded.

Sweep priors (harness.yaml, verbatim from 13_diploshic): s = 10**U(-4,-2); soft f0 recipe-default.
"""
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
HARNESS = f"{D}/scripts/harness"
SLIM_PY = "/home/ssmall/miniforge3/envs/slim_sim/bin/python"
VCF_MODEL = f"{D}/results/vcfretrain_fullsfs/illexModel_vcf"
CHROM = "43"
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1", CUDA_VISIBLE_DEVICES="",
           # cap TensorFlow (diploSHIC predict) threads -> avoid a thread storm across concurrent workers
           TF_NUM_INTRAOP_THREADS="1", TF_NUM_INTEROP_THREADS="1",
           TF_CPP_MIN_LOG_LEVEL="3")
LOCI = list(range(2_000_000, 70_000_000, 2_500_000))   # 28 diverse chr43 starts


def run(cmd, log):
    with open(log, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=ENV)
    return r.returncode


def _trees_path(job, segdir, L):
    if job["klass"] == "neutral":
        return os.path.join(segdir, f"msp_variable_{CHROM}_{job['start']}_s{job['seed']}.recap.trees")
    return os.path.join(segdir, f"{job['klass']}_{CHROM}_{job['start']}_s{job['seed']}.recap.trees")


def do_gen(job, outdir, L, n_dip):
    """Phase 1: forward/coalescent sim -> recapitated tree. CPU-bound; safe at high conc."""
    segdir = os.path.join(outdir, "segs", job["tag"])
    os.makedirs(segdir, exist_ok=True)
    trees = _trees_path(job, segdir, L)
    t0 = time.time()
    if job["klass"] == "neutral":
        gcmd = [SLIM_PY, f"{HARNESS}/msprime_neutral_seg.py", "--chrom", CHROM,
                "--start", str(job["start"]), "--L", str(L), "--seed", str(job["seed"]),
                "--n-dip", str(n_dip), "--recomb", "variable", "--outdir", segdir]
    else:
        gcmd = [SLIM_PY, f"{HARNESS}/calib_run_sim.py", "--klass", job["klass"], "--chrom", CHROM,
                "--start", str(job["start"]), "--L", str(L), "--seed", str(job["seed"]),
                "--n-dip", str(n_dip), "--sel", f"{job['sel']:.6g}", "--outdir", segdir]
    rc = run(gcmd, os.path.join(segdir, "gen.log"))
    ok = (rc == 0 and os.path.exists(trees))
    return dict(job, gen_ok=ok, trees=trees if ok else None, gen_secs=round(time.time() - t0, 1))


def do_scan(job, outdir, L):
    """Phase 2: ascertain -> unphased VCF -> fvecVcf -> predict. Memory-bandwidth-bound; LOW conc."""
    segdir = os.path.join(outdir, "segs", job["tag"])
    trees = job["trees"]
    tag = job["tag"]
    xpos = L // 2
    t0 = time.time()
    scmd = [SLIM_PY, f"{HARNESS}/calib_scan_segment.py", "--trees", trees, "--chrom", CHROM,
            "--start", str(job["start"]), "--L", str(L), "--xpos", str(xpos),
            "--seed", str(job["seed"]), "--path", "vcf", "--model", VCF_MODEL,
            "--outdir", segdir, "--tag", tag]
    rc = run(scmd, os.path.join(segdir, "scan.log"))
    preds = os.path.join(segdir, f"{tag}.vcf.preds")
    ok = os.path.exists(preds)
    if ok:
        try:
            os.remove(trees)     # free the big tree once scanned
        except OSError:
            pass
    return dict(job, status="OK" if ok else "SCAN_FAIL", preds=preds if ok else None,
                xpos=xpos, scan_secs=round(time.time() - t0, 1))


def build_jobs(n_hard, n_soft, n_neutral, rng):
    jobs = []
    def draw_sel():
        return float(10 ** rng.uniform(-4, -2))
    seed = 1000
    for i in range(n_hard):
        seed += 1
        jobs.append(dict(klass="hard", start=LOCI[i % len(LOCI)], seed=seed, sel=draw_sel(),
                         tag=f"hard_{LOCI[i % len(LOCI)]}_s{seed}"))
    for i in range(n_soft):
        seed += 1
        jobs.append(dict(klass="soft", start=LOCI[(i + 7) % len(LOCI)], seed=seed, sel=draw_sel(),
                         tag=f"soft_{LOCI[(i + 7) % len(LOCI)]}_s{seed}"))
    for i in range(n_neutral):
        seed += 1
        st = LOCI[(i * 3) % len(LOCI)]
        jobs.append(dict(klass="neutral", start=st, seed=seed, sel=0.0,
                         tag=f"neutral_{st}_s{seed}"))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-hard", type=int, default=12)
    ap.add_argument("--n-soft", type=int, default=12)
    ap.add_argument("--n-neutral", type=int, default=48)
    ap.add_argument("--L", type=int, default=3_000_000)
    ap.add_argument("--n-dip", type=int, default=350)
    ap.add_argument("--conc-gen", type=int, default=12, help="gen concurrency (CPU-bound)")
    ap.add_argument("--conc-scan", type=int, default=5, help="scan concurrency (memory-bandwidth-bound)")
    ap.add_argument("--seed", type=int, default=20260914)
    ap.add_argument("--outdir", default=f"{D}/calibration/sims")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    import pandas as pd
    rng = np.random.default_rng(a.seed)
    jobs = build_jobs(a.n_hard, a.n_soft, a.n_neutral, rng)
    manifest = os.path.join(a.outdir, "manifest.tsv")
    print(f"[campaign] {len(jobs)} segments (hard={a.n_hard} soft={a.n_soft} neutral={a.n_neutral}) "
          f"L={a.L} n_dip={a.n_dip} conc_gen={a.conc_gen} conc_scan={a.conc_scan} -> {a.outdir}", flush=True)

    # ---- Phase 1: generate all trees (high conc) ----
    genned = []
    done = 0
    with ThreadPoolExecutor(max_workers=a.conc_gen) as ex:
        futs = {ex.submit(do_gen, j, a.outdir, a.L, a.n_dip): j for j in jobs}
        for fut in as_completed(futs):
            r = fut.result(); done += 1
            genned.append(r)
            print(f"[gen] {done}/{len(jobs)} {r['tag']} {'OK' if r['gen_ok'] else 'GEN_FAIL'} {r['gen_secs']}s",
                  flush=True)
    ok_gen = [r for r in genned if r["gen_ok"]]
    print(f"[campaign] gen phase done: {len(ok_gen)}/{len(jobs)} trees", flush=True)

    # ---- Phase 2: scan all trees (low conc, avoids memory-bandwidth wall) ----
    rows = []
    done = 0
    with ThreadPoolExecutor(max_workers=a.conc_scan) as ex:
        futs = {ex.submit(do_scan, j, a.outdir, a.L): j for j in ok_gen}
        for fut in as_completed(futs):
            r = fut.result(); done += 1
            rows.append(r)
            print(f"[scan] {done}/{len(ok_gen)} {r['tag']} {r['status']} {r['scan_secs']}s", flush=True)
            pd.DataFrame(rows).to_csv(manifest, sep="\t", index=False)
    # include gen failures in manifest
    for r in genned:
        if not r["gen_ok"]:
            rows.append(dict(r, status="GEN_FAIL", preds=None))
    df = pd.DataFrame(rows)
    df.to_csv(manifest, sep="\t", index=False)
    print(f"[campaign] DONE. OK={int((df.status=='OK').sum())}/{len(df)} -> {manifest}", flush=True)


if __name__ == "__main__":
    main()
