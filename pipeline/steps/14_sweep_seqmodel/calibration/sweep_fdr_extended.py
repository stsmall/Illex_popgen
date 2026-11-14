#!/usr/bin/env python3
"""Circular-shift nulls for the BGS-robust concordant regions and for Tier-1 calls (review minor 27, M8).

Same null as tier2_concordance_fdr.py: the joint RAiSD/SweepFinder2 stratified-outlier track is shifted circularly
within each chromosome by an independent uniform offset of whole 100-kb windows (10,000 shifts, seed 20261111),
which keeps the number and clustering of the footprint outliers and breaks only their alignment with diploS/HIC.
 (1) High-confidence regions: concordant regions that contain a Tier-1 call and pass the per-region
     background-selection test (bgs_screen_model.test_regions, SWEEP-like = standardised residual <= -0.3 below the
     coding-density expectation), the 11 of Table 4. Both properties belong to the diploS/HIC region and do not depend
     on the footprint track. In each shift, the regions that become concordant by chance and have both are counted;
     FDR = mean null count / observed count. The count without the Tier-1 condition is also reported.
 (2) Tier-1 calls (hmm_decode_masked/tier1_markov_calls.tsv, autosomal): a call is corroborated when one of the
     100-kb windows it spans is a footprint outlier. Observed corroborated calls are compared with the shifted
     nulls; FDR of the corroborated Tier-1 calls = mean null / observed. Tier-1 as a whole still has no error rate.
Usage: sweep_fdr_extended.py OUTDIR [N]
"""
import importlib.util
import sys

import numpy as np
import pandas as pd

A = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
R = f"{A}/results/empirical_scan_fullsfs"
spec = importlib.util.spec_from_file_location("bgs", f"{A}/calibration/bgs_screen_model.py")
bgs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bgs)


def regions(df, gap=300000):
    o = df[df.S_out == 1].sort_values(["chrom", "start"])
    rid, cur, lc, le = [], -1, None, None
    for c, s, e in zip(o.chrom, o.start, o.end):
        if c != lc or s - le > gap:
            cur += 1
        rid.append(cur)
        lc, le = c, e
    return o.index.to_numpy(), np.array(rid)


def main():
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    df = pd.read_csv(f"{R}/outlier_scan_45_masked/windows.tsv", sep="\t", dtype={"chrom": str}) \
        .sort_values(["chrom", "start"]).reset_index(drop=True)
    foot = ((df.SF2_out == 1) | (df.RAiSD_out == 1)).to_numpy()
    idx, rid = regions(df)
    n_reg = rid.max() + 1
    # region extents and BGS verdicts
    rg = pd.DataFrame({"rid": rid, "chrom": df.chrom.to_numpy()[idx], "start": df.start.to_numpy()[idx],
                       "end": df.end.to_numpy()[idx]}).groupby("rid").agg(chrom=("chrom", "first"),
                                                                          start=("start", "min"), end=("end", "max"))
    res = pd.DataFrame(bgs.test_regions([(r.chrom, int(r.start), int(r.end), "") for r in rg.itertuples()]))
    robust = res.verdict.str.startswith("SWEEP-like").to_numpy()

    def hits(track):
        h = np.zeros(n_reg, bool)
        np.logical_or.at(h, rid, track[idx])
        return h
    h0 = hits(foot)
    print(f"diploS/HIC regions {n_reg}, BGS-robust {robust.sum()}; observed concordant {h0.sum()}, "
          f"concordant and BGS-robust {(h0 & robust).sum()}")
    # Tier-1 calls -> window rows
    t1 = pd.read_csv(f"{R}/hmm_decode_masked/tier1_markov_calls.tsv", sep="\t", dtype={"chrom": str})
    t1 = t1[~t1.chrom.isin(["Z"])]
    key = {(c, s): i for i, (c, s) in enumerate(zip(df.chrom, df.start))}
    call_rows = []
    for c, s, e in zip(t1.chrom, t1.start, t1.end):
        rows = [key[(c, p)] for p in range(int(s), int(e), 100000) if (c, p) in key]
        call_rows.append(np.array(rows, int))
    have = np.array([len(r) > 0 for r in call_rows])
    reg_of_row = np.full(len(df), -1)
    reg_of_row[idx] = rid
    has_t1 = np.zeros(n_reg, bool)
    for r in call_rows:
        for row in r:
            if reg_of_row[row] >= 0:
                has_t1[reg_of_row[row]] = True
    hc = robust & has_t1
    print(f"regions with a Tier-1 call {has_t1.sum()}; concordant, Tier-1 and BGS-robust (high-confidence) {(h0 & hc).sum()}")
    def t1_hits(track):
        return np.array([track[r].any() if len(r) else False for r in call_rows])
    o1 = int(t1_hits(foot).sum())
    print(f"Tier-1 autosomal calls {len(t1)} ({have.sum()} mapped to scored windows), corroborated {o1}")
    chrom_rows = [g.index.to_numpy() for _, g in df.groupby("chrom", sort=False)]
    rng = np.random.default_rng(20261111)
    nc, ncr, nhc, nt1 = np.empty(n, int), np.empty(n, int), np.empty(n, int), np.empty(n, int)
    for i in range(n):
        t = np.empty_like(foot)
        for rows in chrom_rows:
            t[rows] = np.roll(foot[rows], rng.integers(len(rows)))
        h = hits(t)
        nc[i], ncr[i], nhc[i], nt1[i] = h.sum(), (h & robust).sum(), (h & hc).sum(), t1_hits(t).sum()
    obs_r = int((h0 & robust).sum())
    obs_hc = int((h0 & hc).sum())
    s = pd.DataFrame([{
        "concordant_obs": int(h0.sum()), "concordant_null_mean": nc.mean(), "concordant_fdr": nc.mean() / h0.sum(),
        "robust_obs": obs_r, "robust_null_mean": ncr.mean(), "robust_null_q975": np.quantile(ncr, 0.975),
        "robust_null_max": ncr.max(), "robust_fdr": ncr.mean() / obs_r, "robust_p": (ncr >= obs_r).mean(),
        "highconf_obs": obs_hc, "highconf_null_mean": nhc.mean(), "highconf_null_q975": np.quantile(nhc, 0.975),
        "highconf_null_max": nhc.max(), "highconf_fdr": nhc.mean() / obs_hc, "highconf_p": (nhc >= obs_hc).mean(),
        "tier1_calls": len(t1), "tier1_corroborated_obs": o1, "tier1_null_mean": nt1.mean(),
        "tier1_null_q975": np.quantile(nt1, 0.975), "tier1_corroborated_fdr": nt1.mean() / o1,
        "tier1_p": (nt1 >= o1).mean(), "n_shifts": n}])
    s.to_csv(f"{out}/sweep_fdr_extended_summary.tsv", sep="\t", index=False)
    res.assign(concordant=h0, has_tier1=has_t1).to_csv(f"{out}/diploshic_regions_bgs.tsv", sep="\t", index=False)
    print(s.T.to_string(header=False))


if __name__ == "__main__":
    main()
