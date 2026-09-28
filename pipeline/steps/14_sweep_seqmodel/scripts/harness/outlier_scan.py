#!/usr/bin/env python3
"""outlier_scan.py -- OPTION 1 sweep scan: empirical-outlier calibration of the fvecVcf-retrained
diploSHIC scan, BGS/accessibility-stratified, with RAiSD+SF2 adding weight (user plan 2026-08-10).

The fvecVcf model still over-calls the real genome (pervasive sim-vs-real diversity-heterogeneity
gap), so we do NOT trust an absolute sim-null. Instead: sweep score S = p_hard+p_soft per window,
rank it WITHIN strata of (callable_frac quintile x cds_frac{0,>0}) so a window only stands out if
its sweep signal is high relative to windows of similar accessibility + gene density (controls the
accessibility+BGS baseline that inflates S). RAiSD (maxRAiSD_mu) and SF2 (maxSF2_LR), reused from
the old integration (model-independent), are stratum-ranked the same way and add concordance weight.
Adjacent diploSHIC-outlier windows are merged (recomb-aware gap via the ReLERNN map_df, else physical)
into candidate sweep REGIONS, ranked by combined evidence and split hard vs soft.

Usage: outlier_scan.py --preds genome.preds --integ integrated.tsv [--mapdf mapdf.tsv]
       [--q 0.99] [--min-callable 0.25] [--merge-gap-cm 0.05] --out OUTDIR
"""
import argparse, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))  # for seqmodel if needed later
D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
SCHEMA = f"{D}/config/probability_schema.yaml"


def load_preds(path):
    sys.path.insert(0, f"{D}/scripts")
    from seqmodel.std_diploshic import standardize
    df = standardize(path, SCHEMA)             # chrom,start,end,p_neutral,p_linked_soft,p_linked_hard,p_soft,p_hard
    df["S"] = df["p_hard"] + df["p_soft"]      # sweep score (center)
    df["chrom"] = df["chrom"].astype(str)
    return df


def strat_pct(df, val, strata):
    """within-stratum percentile rank (0-1) of column `val`, grouped by `strata` cols."""
    return df.groupby(strata)[val].rank(pct=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preds", required=True); ap.add_argument("--integ", required=True)
    ap.add_argument("--mapdf", default=f"{D}/config/mapdf.tsv")
    ap.add_argument("--q", type=float, default=0.99); ap.add_argument("--min-callable", type=float, default=0.25)
    ap.add_argument("--merge-gap-bp", type=int, default=300000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mask", default=None, help="BED of regions to exclude before stratified ranking (e.g. collapsed multi-copy blocks)")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

    p = load_preds(a.preds)
    integ = pd.read_csv(a.integ, sep="\t")
    integ["chrom"] = integ["chrom"].astype(str)
    # join by 100kb window index (preds classifiedWinStart is 1-based; integ start is 0-based)
    p["wkey"] = (p["start"] - 1) // 100000
    integ["wkey"] = integ["start"] // 100000
    keep = ["chrom", "wkey", "maxSF2_LR", "maxRAiSD_mu", "callable_frac", "cds_frac"]
    df = p.merge(integ[keep], on=["chrom", "wkey"], how="left")
    miss = df["callable_frac"].isna().mean()
    if miss > 0.05:
        sys.stderr.write(f"[warn] {100*miss:.0f}% windows unmatched to integ (check coord join)\n")
    for c in ["maxSF2_LR", "maxRAiSD_mu", "callable_frac", "cds_frac"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df[df["callable_frac"] >= a.min_callable].reset_index(drop=True)
    if a.mask:
        m = pd.read_csv(a.mask, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "s", "e"]); m["chrom"] = m["chrom"].astype(str)
        drop = np.zeros(len(df), bool)
        for r in m.itertuples():
            drop |= (df["chrom"].values == r.chrom) & (df["start"].values < r.e) & (df["end"].values > r.s)
        print(f"[mask] excluded {drop.sum()} windows in {len(m)} masked regions"); df = df[~drop].reset_index(drop=True)

    # strata: callable quintile x cds {0,>0}
    df["call_q"] = pd.qcut(df["callable_frac"].rank(method="first"), 5, labels=False)
    df["cds_bin"] = (df["cds_frac"].fillna(0) > 0).astype(int)
    strata = ["call_q", "cds_bin"]
    df["S_pct"] = strat_pct(df, "S", strata)
    df["SF2_pct"] = strat_pct(df, "maxSF2_LR", strata)
    df["RAiSD_pct"] = strat_pct(df, "maxRAiSD_mu", strata)
    df["S_out"] = (df["S_pct"] >= a.q).astype(int)
    df["SF2_out"] = (df["SF2_pct"] >= a.q).astype(int)
    df["RAiSD_out"] = (df["RAiSD_pct"] >= a.q).astype(int)
    df["n_methods"] = df[["S_out", "SF2_out", "RAiSD_out"]].sum(1)
    # hard vs soft label for diploSHIC outliers
    df["sweep_type"] = np.where(df["p_hard"] >= df["p_soft"], "hard", "soft")
    df.to_csv(f"{a.out}/windows.tsv", sep="\t", index=False)

    # ---- merge adjacent diploSHIC-S outlier windows -> regions ----
    o = df[df["S_out"] == 1].sort_values(["chrom", "start"]).reset_index(drop=True)
    regions = []
    if len(o):
        cur = None
        for _, w in o.iterrows():
            if cur and w["chrom"] == cur["chrom"] and w["start"] - cur["end"] <= a.merge_gap_bp:
                cur["end"] = w["end"]; cur["rows"].append(w)
            else:
                if cur: regions.append(cur)
                cur = {"chrom": w["chrom"], "start": w["start"], "end": w["end"], "rows": [w]}
        if cur: regions.append(cur)
    rrows = []
    for r in regions:
        R = pd.DataFrame(r["rows"])
        rrows.append(dict(chrom=r["chrom"], start=r["start"], end=r["end"],
                          n_win=len(R), maxS=R["S"].max(), maxS_pct=R["S_pct"].max(),
                          sweep_type=R.loc[R["S"].idxmax(), "sweep_type"],
                          sf2_support=int((R["SF2_out"] == 1).any()), raisd_support=int((R["RAiSD_out"] == 1).any()),
                          max_nmethods=int(R["n_methods"].max()),
                          mean_cds=R["cds_frac"].mean(), mean_callable=R["callable_frac"].mean()))
    reg = pd.DataFrame(rrows)
    if len(reg):
        reg["evidence"] = reg["maxS_pct"] + reg["sf2_support"] + reg["raisd_support"]  # diploSHIC + concordance
        reg = reg.sort_values(["max_nmethods", "evidence"], ascending=False).reset_index(drop=True)
    reg.to_csv(f"{a.out}/regions.tsv", sep="\t", index=False)

    # ---- summary ----
    n = len(df)
    print(f"[outlier] {n} callable windows (>= {a.min_callable} callable)")
    print(f"[outlier] diploSHIC-S outliers (top {100*(1-a.q):.0f}% within stratum): {df['S_out'].sum()} "
          f"({100*df['S_out'].mean():.2f}% of windows)  hard={int((df.S_out&(df.sweep_type=='hard')).sum())} "
          f"soft={int((df.S_out&(df.sweep_type=='soft')).sum())}")
    print(f"[outlier] candidate REGIONS: {len(reg)}")
    if len(reg):
        conc = reg[reg["max_nmethods"] >= 2]
        print(f"[outlier]   with >=2-method concordance (diploSHIC+RAiSD/SF2): {len(conc)}")
        print(f"[outlier]   with 3-method: {int((reg['max_nmethods']==3).sum())}")
        print("\nTop 15 candidate regions (by concordance, then evidence):")
        cols = ["chrom", "start", "end", "n_win", "sweep_type", "maxS", "maxS_pct",
                "sf2_support", "raisd_support", "max_nmethods", "mean_cds"]
        print(reg[cols].head(15).to_string(index=False))
    print(f"\n[outlier] wrote {a.out}/windows.tsv + regions.tsv")


if __name__ == "__main__":
    main()
