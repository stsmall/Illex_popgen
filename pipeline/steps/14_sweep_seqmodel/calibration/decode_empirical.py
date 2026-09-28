#!/usr/bin/env python3
"""decode_empirical.py -- PLAN_02 Phase E: decode the empirical per-window diploSHIC
emissions across all 45 chroms with the recomb-aware HMM (and HSMM) using the calibrated
params in calibration/params.json, and emit Markov-corrected sweep-call REGIONS.

Input : results/empirical_scan_fullsfs/outlier_scan_45/windows.tsv  (already-standardized
        per-100kb-window class probs p_neutral,p_linked_hard,p_linked_soft,p_hard,p_soft
        for all 45 chroms incl chr2 & chr42).
Map   : config/mapdf.tsv (sex-avg ReLERNN). chr2 & chr42 absent -> constant-rate fallback
        at the genome length-weighted mean rate (documented in params.json).
Params: calibration/params.json (emission transform, HMM transition, HSMM duration, null
        thresholds). Null maxima arrays (for empirical p) read from calibration/null_*.npy.

Output: <out>/regions.tsv        (chrom,start,end,score,state,evidence + p-value + pass flag)
        <out>/state_path.tsv.gz   (per-window Viterbi state, HMM)
"""
import argparse, json, os, sys
import numpy as np, pandas as pd

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
sys.path.insert(0, f"{D}/scripts")
from seqmodel.emissions import collapse, emission_logprob
from seqmodel.mapcoord import split_blocks, deltas, genetic_positions
from seqmodel.decode import viterbi_hmm, viterbi_hsmm
from seqmodel.regions import extract_regions
from seqmodel.calibrate_null import empirical_p

PROB = ["p_neutral", "p_linked_hard", "p_linked_soft", "p_hard", "p_soft"]


def build_map(mapdf_path, windows, fallback_rate):
    """Return a map_df (chrom,start,end,rate). For chroms absent from mapdf (chr2/chr42),
    synthesize a single constant-rate interval spanning that chrom's windows."""
    m = pd.read_csv(mapdf_path, sep="\t")
    m["chrom"] = m["chrom"].astype(str)
    have = set(m["chrom"])
    extra = []
    for c, g in windows.groupby(windows["chrom"].astype(str)):
        if c not in have:
            lo = int(g["start"].min()) - 1
            hi = int(g["end"].max()) + 1
            extra.append(dict(chrom=c, start=max(lo, 0), end=hi, rate=fallback_rate))
    if extra:
        m = pd.concat([m, pd.DataFrame(extra)], ignore_index=True)
    return m[["chrom", "start", "end", "rate"]]


def decode_chrom(wc, map_df, chrom, P, model, max_phys_gap):
    """Decode one chromosome. Returns (path list len n, regions_df, collapsed, emis_log)."""
    emis_log = emission_logprob(collapse(wc), transform=P["emission"]["transform"],
                                priors=tuple(P["emission"]["priors"]),
                                temperature=P["emission"].get("temperature", 1.0))
    col = collapse(wc)
    blocks = split_blocks(wc, map_df, chrom, max_phys_gap)
    path = ["GAP"] * len(wc)
    gpos_full = np.full(len(wc), np.nan)
    for (a, b) in blocks:
        wb = wc.iloc[a:b + 1]
        if model == "hmm":
            dg, dx = deltas(wb, map_df, chrom) if b > a else (np.array([]), np.array([]))
            p, _ = viterbi_hmm(emis_log[a:b + 1], dg, dx, P["hmm"])
        else:  # hsmm
            gp = genetic_positions(wb, map_df, chrom)
            gpos_full[a:b + 1] = gp
            p, _ = viterbi_hsmm(emis_log[a:b + 1], gp - gp[0], P["hsmm"])
        path[a:b + 1] = p
    # region scores need dg/dx (hmm) or gpos (hsmm) over the WHOLE chrom index space
    if model == "hmm":
        # per-consecutive-window dg/dx across the chrom (regions are within blocks, contiguous)
        gfull = genetic_positions(wc, map_df, chrom)
        dg_full = np.diff(gfull)
        dx_full = np.diff(((wc["start"].to_numpy(np.int64) + wc["end"].to_numpy(np.int64)) // 2).astype(float))
        R = extract_regions(wc, path, emis_log, col, params=P["hmm"], model="hmm",
                            dg=dg_full, dx=dx_full)
    else:
        gfull = genetic_positions(wc, map_df, chrom)
        R = extract_regions(wc, path, emis_log, col, params=P["hsmm"], model="hsmm",
                            gpos=gfull)
    return path, R, col, emis_log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", default=f"{D}/results/empirical_scan_fullsfs/outlier_scan_45/windows.tsv")
    ap.add_argument("--mapdf", default=f"{D}/config/mapdf.tsv")
    ap.add_argument("--params", default=f"{D}/calibration/params.json")
    ap.add_argument("--out", default=f"{D}/results/empirical_scan_fullsfs/hmm_decode")
    ap.add_argument("--model", choices=["hmm", "hsmm", "both"], default="both")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    P = json.load(open(a.params))
    fallback = P["decode"]["fallback_rate"]
    max_phys_gap = P["decode"]["max_phys_gap"]

    win = pd.read_csv(a.windows, sep="\t")
    win["chrom"] = win["chrom"].astype(str)
    win = win.sort_values(["chrom", "start"]).reset_index(drop=True)
    map_df = build_map(a.mapdf, win, fallback)

    def chrom_key(c):
        return (0, int(c)) if c.isdigit() else (1, c)

    models = ["hmm", "hsmm"] if a.model == "both" else [a.model]
    for model in models:
        null_path = os.path.join(os.path.dirname(a.params), f"null_maxima_{model}.npy")
        null_max = np.load(null_path) if os.path.exists(null_path) else None
        thr = P["null"].get(f"threshold_{model}", None)
        all_regions = []
        all_paths = []
        for chrom in sorted(win["chrom"].unique(), key=chrom_key):
            wc = win[win["chrom"] == chrom].reset_index(drop=True)
            path, R, col, emis = decode_chrom(wc, map_df, chrom, P, model, max_phys_gap)
            pdf = pd.DataFrame({"chrom": chrom, "start": wc["start"], "end": wc["end"], "state": path})
            all_paths.append(pdf)
            if len(R):
                R = R.copy()
                R["model"] = model
                R["chrom_fallback"] = chrom in ("2", "42")   # constant-rate map (absent from ReLERNN mapdf)
                all_regions.append(R)
        reg = pd.concat(all_regions, ignore_index=True) if all_regions else pd.DataFrame()
        if len(reg):
            if null_max is not None:
                reg["p_value"] = [empirical_p(s, null_max) for s in reg["score"]]
            if thr is not None:
                reg["pass_threshold"] = reg["score"] >= thr
            # evidence label: hard vs soft from evidence sums
            reg["evidence"] = np.where(reg["hard_evidence"] >= reg["soft_evidence"], "hard", "soft")
            reg = reg.sort_values("score", ascending=False).reset_index(drop=True)
        out_reg = os.path.join(a.out, f"regions_{model}.tsv")
        reg.to_csv(out_reg, sep="\t", index=False)
        paths = pd.concat(all_paths, ignore_index=True)
        paths.to_csv(os.path.join(a.out, f"state_path_{model}.tsv.gz"), sep="\t", index=False,
                     compression="gzip")
        npass = int(reg["pass_threshold"].sum()) if (len(reg) and "pass_threshold" in reg) else 0
        print(f"[decode:{model}] {len(reg)} candidate regions, {npass} pass threshold "
              f"(thr={thr}); wrote {out_reg}")
        # primary deliverable path for HMM
        if model == "hmm":
            reg.to_csv(os.path.join(a.out, "regions.tsv"), sep="\t", index=False)
            paths.to_csv(os.path.join(a.out, "state_path.tsv.gz"), sep="\t", index=False,
                         compression="gzip")


if __name__ == "__main__":
    main()
