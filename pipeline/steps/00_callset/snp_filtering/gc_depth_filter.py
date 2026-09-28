#!/usr/bin/env python3
"""
Compute GC-normalized depth mask.
Flags genomic windows where depth exceeds GC-matched expectation by > N SD.

Usage:
    python3 gc_depth_filter.py \
        --depth sample.regions.bed.gz \
        --gc windows_gc.bed \
        --output high_depth_windows.bed \
        --sd-threshold 2.0
"""

import argparse
import pandas as pd
import numpy as np
import statsmodels.api as sm


def main():
    parser = argparse.ArgumentParser(
        description="Flag genomic windows with GC-normalized excess depth"
    )
    parser.add_argument("--depth", required=True,
                        help="mosdepth regions BED.gz (chrom start end depth)")
    parser.add_argument("--gc", required=True,
                        help="bedtools nuc output: chrom start end gc_frac")
    parser.add_argument("--output", required=True,
                        help="Output BED of high-depth windows")
    parser.add_argument("--sd-threshold", type=float, default=2.0,
                        help="SD above GC-matched expectation to flag (default: 2.0)")
    parser.add_argument("--lowess-frac", type=float, default=0.3,
                        help="Fraction of data for LOWESS bandwidth (default: 0.3)")
    parser.add_argument("--min-depth", type=float, default=0.0,
                        help="Exclude windows with depth <= this before fitting (default: 0)")
    args = parser.parse_args()

    # Load depth
    depth = pd.read_csv(
        args.depth, sep="\t", compression="infer",
        names=["chrom", "start", "end", "depth"]
    )

    # Load GC content
    gc = pd.read_csv(
        args.gc, sep="\t",
        names=["chrom", "start", "end", "gc_frac"]
    )

    df = depth.merge(gc, on=["chrom", "start", "end"])

    # Remove zero/low-depth windows before fitting to avoid distorting the curve
    df = df[df["depth"] > args.min_depth].copy()

    if df.empty:
        raise ValueError("No windows remain after depth filtering — check inputs")

    print(f"  Windows for fitting: {len(df)}")
    print(f"  GC range: {df['gc_frac'].min():.3f} – {df['gc_frac'].max():.3f}")
    print(f"  Depth range: {df['depth'].min():.1f} – {df['depth'].max():.1f}")

    # LOWESS fit: depth ~ GC content
    lowess = sm.nonparametric.lowess(
        df["depth"],
        df["gc_frac"],
        frac=args.lowess_frac,
        return_sorted=True
    )

    # Interpolate expected depth for each window from LOWESS curve
    df["expected_depth"] = np.interp(
        df["gc_frac"],
        lowess[:, 0],
        lowess[:, 1]
    )
    df["depth_residual"] = df["depth"] - df["expected_depth"]

    residual_mean = df["depth_residual"].mean()
    residual_sd   = df["depth_residual"].std()
    threshold     = residual_mean + args.sd_threshold * residual_sd

    flagged = df[df["depth_residual"] > threshold][["chrom", "start", "end"]]
    flagged.to_csv(args.output, sep="\t", index=False, header=False)

    n_total   = len(df)
    n_flagged = len(flagged)
    print(f"  Depth residual mean: {residual_mean:.2f}")
    print(f"  Depth residual SD:   {residual_sd:.2f}")
    print(f"  Flagging threshold:  {threshold:.2f} ({args.sd_threshold} SD above mean)")
    print(f"  Windows flagged:     {n_flagged} / {n_total} ({100*n_flagged/n_total:.1f}%)")
    print(f"  Output:              {args.output}")


if __name__ == "__main__":
    main()
