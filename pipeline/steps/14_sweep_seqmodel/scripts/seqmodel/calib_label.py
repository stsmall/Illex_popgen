#!/usr/bin/env python3
"""calib_label.py -- PLAN_02 Phase C: label tiled calibration windows N/LL/C/LR by GENETIC
distance from the known sweep centre XPOS, and turn a per-segment slice of the genome map
into a decoder map_df in the segment's local ("seg") coordinates.

A window is Center (C) if its midpoint is within d_C Morgans of XPOS; Linked (LL left / LR
right) if within (d_C, d_L]; Neutral (N) beyond d_L. Distances use the SAME recombination map
that generated the sim, so the label reflects genetic — not physical — proximity (the whole
point of the recomb-aware decoder). d_C / d_L are tuned in Phase C; this module is parameterized.
"""
import numpy as np
import pandas as pd


def seg_mapdf(genome_mapdf, chrom, start, L, seg_name="seg"):
    """Slice the genome map_df to [start, start+L) on `chrom` and shift to local coords
    (0..L), renaming chrom -> seg_name. Matches the synthetic 'seg' contig the calib scan
    tiles, so mapcoord covers every tiled-window midpoint the sim actually produced."""
    g = genome_mapdf[genome_mapdf["chrom"].astype(str) == str(chrom)].copy()
    g = g[(g["end"] > start) & (g["start"] < start + L)]
    g["start"] = np.clip(g["start"] - start, 0, L)
    g["end"] = np.clip(g["end"] - start, 0, L)
    g["chrom"] = seg_name
    return g.sort_values("start").reset_index(drop=True)[["chrom", "start", "end", "rate"]]


def _genetic_at(map_df, bp):
    """Cumulative genetic position (Morgans) at physical bp using contiguous map_df intervals."""
    s = map_df["start"].to_numpy(np.int64); e = map_df["end"].to_numpy(np.int64)
    r = map_df["rate"].to_numpy(float)
    cum = np.concatenate(([0.0], np.cumsum((e[:-1] - s[:-1]) * r[:-1])))
    idx = int(np.searchsorted(s, bp, side="right") - 1)
    idx = min(max(idx, 0), len(s) - 1)
    return cum[idx] + (bp - s[idx]) * r[idx]


def label_windows(preds_df, map_df, xpos, d_C, d_L, seg_name="seg"):
    """Add gdist_cM (signed cM from xpos; +=right of sweep) and region (N/LL/C/LR) columns.

    preds_df: standardized (chrom,start,end,p_*). map_df: seg-coord contiguous intervals.
    d_C, d_L: Morgans. Windows whose midpoint the map doesn't cover -> region 'GAP'."""
    df = preds_df.copy()
    mids = ((df["start"].to_numpy(np.int64) + df["end"].to_numpy(np.int64)) // 2)
    gx = _genetic_at(map_df, int(xpos))
    smin, emax = int(map_df["start"].min()), int(map_df["end"].max())
    reg, gd = [], []
    for m in mids:
        if m < smin or m >= emax:
            reg.append("GAP"); gd.append(np.nan); continue
        gm = _genetic_at(map_df, int(m))
        signed = gm - gx
        dist = abs(signed)
        if dist <= d_C:
            reg.append("C")
        elif dist <= d_L:
            reg.append("LL" if signed < 0 else "LR")
        else:
            reg.append("N")
        gd.append(signed * 100.0)                       # cM
    df["gdist_cM"] = gd
    df["region"] = reg
    return df


def summarize(labeled_df):
    """One row per window: coord, predClass (argmax), region label, gdist — for eyeballing."""
    PROB = ["p_neutral", "p_linked_soft", "p_linked_hard", "p_soft", "p_hard"]
    short = {"p_neutral": "neut", "p_linked_soft": "lSoft", "p_linked_hard": "lHard",
             "p_soft": "soft", "p_hard": "hard"}
    P = labeled_df[PROB].to_numpy(float)
    pred = [short[PROB[i]] for i in P.argmax(1)]
    out = labeled_df[["start", "end", "region", "gdist_cM"]].copy()
    out["pred"] = pred
    out["pC"] = (labeled_df["p_hard"] + labeled_df["p_soft"]).round(3)
    out["pL"] = (labeled_df["p_linked_hard"] + labeled_df["p_linked_soft"]).round(3)
    out["pN"] = labeled_df["p_neutral"].round(3)
    return out
