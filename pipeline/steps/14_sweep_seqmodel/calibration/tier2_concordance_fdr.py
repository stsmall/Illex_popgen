#!/usr/bin/env python3
"""tier2_concordance_fdr.py -- false-discovery rate of the cross-method-concordant (Tier-2) sweep regions.

Tier-2 regions are diploS/HIC stratified-outlier regions (S_out windows merged within 300 kb by outlier_scan.py)
that contain at least one window that is also a stratified outlier for RAiSD or SweepFinder2 (n_methods >= 2).
The null keeps the diploS/HIC outlier regions and the joint RAiSD/SweepFinder2 outlier track as they are, and
breaks only their alignment: the footprint track is shifted circularly within each chromosome by an independent
uniform offset (whole 100-kb windows), which preserves the number, clustering and RAiSD-SF2 correlation of the
footprint outliers and the per-chromosome outlier density. The number of regions concordant by chance is counted
in each of --n shifts. The FDR of the concordant set is E[null concordant] / observed concordant.

Usage: tier2_concordance_fdr.py --windows .../outlier_scan_45_masked/windows.tsv --out OUTDIR [--n 10000]
Writes OUTDIR/tier2_concordance_fdr.tsv (one row per shift) and OUTDIR/tier2_concordance_fdr_summary.tsv.
"""
import argparse
import os

import numpy as np
import pandas as pd


def regions(df, gap=300000):
    """Region id per diploS/HIC outlier window, merged as in outlier_scan.py (same chromosome, gap <= 300 kb)."""
    o = df[df.S_out == 1].sort_values(["chrom", "start"])
    rid, cur, last_chrom, last_end = [], -1, None, None
    for c, s, e in zip(o.chrom, o.start, o.end):
        if c != last_chrom or s - last_end > gap:
            cur += 1
        rid.append(cur)
        last_chrom, last_end = c, e
    return o.index.to_numpy(), np.array(rid)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261111)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    df = pd.read_csv(a.windows, sep="\t", dtype={"chrom": str}).sort_values(["chrom", "start"]).reset_index(drop=True)
    foot = ((df.SF2_out == 1) | (df.RAiSD_out == 1)).to_numpy()
    idx, rid = regions(df)
    n_reg = rid.max() + 1

    def n_conc(track):
        hit = np.zeros(n_reg, bool)
        np.logical_or.at(hit, rid, track[idx])
        return int(hit.sum())

    obs = n_conc(foot)
    # windows are contiguous 100-kb bins within each chromosome after the callable filter; shift within chromosome
    chrom_rows = [g.index.to_numpy() for _, g in df.groupby("chrom", sort=False)]
    rng = np.random.default_rng(a.seed)
    null = np.empty(a.n, int)
    for i in range(a.n):
        t = np.empty_like(foot)
        for rows in chrom_rows:
            t[rows] = np.roll(foot[rows], rng.integers(len(rows)))
        null[i] = n_conc(t)
    pd.DataFrame({"shift": np.arange(a.n), "n_concordant": null}).to_csv(
        f"{a.out}/tier2_concordance_fdr.tsv", sep="\t", index=False)
    summ = pd.DataFrame([{
        "windows": len(df), "diploshic_outlier_windows": int(df.S_out.sum()), "diploshic_regions": int(n_reg),
        "footprint_outlier_windows": int(foot.sum()), "observed_concordant": obs,
        "null_mean": float(null.mean()), "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)), "null_max": int(null.max()),
        "p_null_ge_obs": float((null >= obs).mean()), "fdr": float(null.mean() / obs) if obs else np.nan,
        "n_shifts": a.n, "seed": a.seed}])
    summ.to_csv(f"{a.out}/tier2_concordance_fdr_summary.tsv", sep="\t", index=False)
    print(summ.T.to_string(header=False))


if __name__ == "__main__":
    main()
