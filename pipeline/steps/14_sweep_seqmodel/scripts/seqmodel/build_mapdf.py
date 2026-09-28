#!/usr/bin/env python3
"""build_mapdf.py -- genome-wide sex-averaged recombination map_df for the HMM decoder.

Reads the ReLERNN male + female BSCORRECTED autosomal maps and produces a contiguous,
non-overlapping, gap-filled per-chromosome map (chrom, start, end, rate) that
`seqmodel.mapcoord` consumes (rate = Morgans/bp = ReLERNN recombRate). Sex-averaged =
mean of male+female where both cover a base, else whichever covers; internal gaps (windows
ReLERNN didn't estimate) are filled with the chromosome's length-weighted mean rate so every
scan-window midpoint inside the ReLERNN span is covered (windows OUTSIDE the span fall in no
interval -> the decoder's split_blocks drops them to GAP, which is correct: no estimate there).

This mirrors relernn_to_slim_map.py's sex-average + gap-fill, but at INTERVAL resolution
(no per-bp array) so it scales to the whole genome. chrom col is a bytes-repr ("b'43'")->stripped.

Usage: build_mapdf.py --male <M.BSCORRECTED.txt> --female <F.BSCORRECTED.txt> --out mapdf.tsv
"""
import argparse
import numpy as np
import pandas as pd


def load_map(path):
    df = pd.read_csv(path, sep="\t")
    df["chrom"] = df["chrom"].astype(str).str.replace(r"^b'|'$", "", regex=True)
    return df[["chrom", "start", "end", "recombRate"]].rename(columns={"recombRate": "rate"})


def _sex_rate_at(g, mids):
    """rate at each midpoint from one sex's intervals g (sorted by start); NaN if uncovered."""
    if len(g) == 0:
        return np.full(len(mids), np.nan)
    s = g.start.to_numpy(np.int64); e = g.end.to_numpy(np.int64); r = g.rate.to_numpy(float)
    idx = np.searchsorted(s, mids, side="right") - 1
    ok = (idx >= 0) & (mids < e[np.clip(idx, 0, len(e) - 1)])
    return np.where(ok, r[np.clip(idx, 0, len(r) - 1)], np.nan)


def build_chrom(m, f):
    """Contiguous gap-filled sex-averaged intervals for one chrom from male df m, female df f."""
    bps = np.unique(np.concatenate([
        m.start.to_numpy(np.int64), m.end.to_numpy(np.int64),
        f.start.to_numpy(np.int64), f.end.to_numpy(np.int64)]))
    if len(bps) < 2:
        return []
    lo, hi = bps[:-1], bps[1:]
    mids = (lo + hi) // 2
    rm = _sex_rate_at(m.sort_values("start"), mids)
    rf = _sex_rate_at(f.sort_values("start"), mids)
    rate = np.nanmean(np.vstack([rm, rf]), axis=0)          # sex-average; NaN only if neither covers
    lengths = (hi - lo).astype(float)
    valid = ~np.isnan(rate)
    fill = float(np.average(rate[valid], weights=lengths[valid])) if valid.any() else 2.1e-9
    rate = np.where(np.isnan(rate), fill, rate)
    # collapse identical-rate consecutive elementary intervals
    rows = []
    cs, ce, cr = lo[0], hi[0], rate[0]
    for i in range(1, len(lo)):
        if rate[i] == cr and lo[i] == ce:
            ce = hi[i]
        else:
            rows.append((cs, ce, cr)); cs, ce, cr = lo[i], hi[i], rate[i]
    rows.append((cs, ce, cr))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--male", required=True)
    ap.add_argument("--female", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    M, F = load_map(a.male), load_map(a.female)
    chroms = sorted(set(M.chrom) | set(F.chrom), key=lambda c: (len(c), c))
    out = []
    for c in chroms:
        rows = build_chrom(M[M.chrom == c], F[F.chrom == c])
        for (s, e, r) in rows:
            out.append((c, s, e, r))
    df = pd.DataFrame(out, columns=["chrom", "start", "end", "rate"])
    df.to_csv(a.out, sep="\t", index=False)
    # summary
    for c in chroms:
        g = df[df.chrom == c]
        span = g.end.max() - g.start.min()
        wmean = np.average(g.rate, weights=(g.end - g.start))
        print(f"chr{c}: {len(g)} intervals, span {g.start.min()}-{g.end.max()} "
              f"({span/1e6:.1f}Mb), mean {wmean*1e8:.3f} cM/Mb")


if __name__ == "__main__":
    main()
