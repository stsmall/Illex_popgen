#!/usr/bin/env python3
"""scan_existing.py -- scan the already-generated calibration trees (calibration/sims/segs/*/
*.recap.trees) through the empirical fvecVcf path (calib_scan_segment.py + illexModel_vcf), at
HIGH concurrency. The ascertainment is compute-bound and parallelizes well (100% CPU/proc even at
conc=12), so conc~15 finishes ~34 segments in ~1.5 h instead of ~5 h at conc=5. Writes manifest.tsv
in the exact shape fit_params.py consumes (klass, start, seed, tag, status, preds, xpos).
"""
import argparse, glob, os, re, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
HARNESS = f"{D}/scripts/harness"
SLIM_PY = "/home/ssmall/miniforge3/envs/slim_sim/bin/python"
VCF_MODEL = f"{D}/results/vcfretrain_fullsfs/illexModel_vcf"
CHROM = "43"
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1", CUDA_VISIBLE_DEVICES="", TF_NUM_INTRAOP_THREADS="1",
           TF_NUM_INTEROP_THREADS="1", TF_CPP_MIN_LOG_LEVEL="3")


def parse_tree(path):
    """(klass, start, seed, tag) from a *.recap.trees filename + its segdir."""
    fn = os.path.basename(path)
    segdir = os.path.dirname(path)
    tag = os.path.basename(segdir)
    if fn.startswith("msp_variable_"):
        m = re.match(r"msp_variable_43_(\d+)_s(\d+)\.recap\.trees", fn)
        klass, start, seed = "neutral", int(m.group(1)), int(m.group(2))
    else:
        m = re.match(r"(hard|soft)_43_(\d+)_s(\d+)\.recap\.trees", fn)
        klass, start, seed = m.group(1), int(m.group(2)), int(m.group(3))
    return klass, start, seed, tag, segdir


def do_scan(path, L):
    klass, start, seed, tag, segdir = parse_tree(path)
    xpos = L // 2
    t0 = time.time()
    scmd = [SLIM_PY, f"{HARNESS}/calib_scan_segment.py", "--trees", path, "--chrom", CHROM,
            "--start", str(start), "--L", str(L), "--xpos", str(xpos), "--seed", str(seed),
            "--path", "vcf", "--model", VCF_MODEL, "--outdir", segdir, "--tag", tag]
    with open(os.path.join(segdir, "scan.log"), "w") as fh:
        rc = subprocess.run(scmd, stdout=fh, stderr=subprocess.STDOUT, env=ENV).returncode
    preds = os.path.join(segdir, f"{tag}.vcf.preds")
    ok = os.path.exists(preds)
    return dict(klass=klass, start=start, seed=seed, tag=tag,
                status="OK" if ok else "SCAN_FAIL", preds=preds if ok else None,
                xpos=xpos, scan_secs=round(time.time() - t0, 1), rc=rc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--segs", default=f"{D}/calibration/sims/segs")
    ap.add_argument("--L", type=int, default=2_500_000)
    ap.add_argument("--conc", type=int, default=15)
    ap.add_argument("--manifest", default=f"{D}/calibration/sims/manifest.tsv")
    a = ap.parse_args()
    trees = sorted(glob.glob(os.path.join(a.segs, "*", "*.recap.trees")))
    print(f"[scan] {len(trees)} trees, conc={a.conc}, L={a.L}", flush=True)
    rows, done = [], 0
    with ThreadPoolExecutor(max_workers=a.conc) as ex:
        futs = {ex.submit(do_scan, t, a.L): t for t in trees}
        for fut in as_completed(futs):
            r = fut.result(); done += 1; rows.append(r)
            print(f"[scan] {done}/{len(trees)} {r['tag']} {r['status']} {r['scan_secs']}s", flush=True)
            pd.DataFrame(rows).to_csv(a.manifest, sep="\t", index=False)
    df = pd.DataFrame(rows)
    df.to_csv(a.manifest, sep="\t", index=False)
    print(f"[scan] DONE OK={int((df.status=='OK').sum())}/{len(df)} "
          f"(hard={int(((df.klass=='hard')&(df.status=='OK')).sum())} "
          f"soft={int(((df.klass=='soft')&(df.status=='OK')).sum())} "
          f"neutral={int(((df.klass=='neutral')&(df.status=='OK')).sum())}) -> {a.manifest}", flush=True)


if __name__ == "__main__":
    main()
