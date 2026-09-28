# scripts/seqmodel/mapcoord.py
import numpy as np
def _check_single_chrom(windows_df, chrom):
    u = [str(x) for x in windows_df["chrom"].unique()]
    if u != [str(chrom)]:
        raise ValueError(f"windows_df must contain exactly chromosome {chrom!r}, found {u}")
def _prep(map_df, chrom):
    g = map_df[map_df.chrom == str(chrom)].sort_values("start")
    s,e,r = g.start.to_numpy(np.int64), g.end.to_numpy(np.int64), g.rate.to_numpy(float)
    if len(s) == 0:
        raise ValueError(f"no map intervals for chromosome {chrom}")
    cum_start = np.concatenate(([0.0], np.cumsum((e[:-1]-s[:-1]) * r[:-1])))  # Morgans at each interval start
    return s,e,r,cum_start
def _mids(w): return ((w.start.to_numpy(np.int64)+w.end.to_numpy(np.int64))//2)
def genetic_positions(windows_df, map_df, chrom):
    _check_single_chrom(windows_df, chrom)
    s,e,r,cum = _prep(map_df, chrom); mids=_mids(windows_df)
    idx = np.searchsorted(s,mids,side="right")-1
    if (idx<0).any() or (mids>=e[np.clip(idx,0,len(e)-1)]).any():
        bad=mids[(idx<0)|(mids>=e[np.clip(idx,0,len(e)-1)])][:5]
        raise ValueError(f"map does not cover midpoint(s) on {chrom}: {list(bad)}")
    return cum[idx] + (mids - s[idx]) * r[idx]
def deltas(windows_df, map_df, chrom):
    g = genetic_positions(windows_df, map_df, chrom); mids=_mids(windows_df).astype(float)
    return np.diff(g), np.diff(mids)
def split_blocks(windows_df, map_df, chrom, max_phys_gap):
    _check_single_chrom(windows_df, chrom)
    mids=_mids(windows_df); s,e,r,_=_prep(map_df,chrom)
    idx=np.searchsorted(s,mids,side="right")-1
    covered=(idx>=0)&(mids<e[np.clip(idx,0,len(e)-1)])
    covered &= r[np.clip(idx,0,len(r)-1)] > 0   # zero-rate interval => break
    gap_after = np.zeros(len(s), dtype=bool)
    gap_after[:-1] = s[1:] != e[:-1]          # True where interval k and k+1 are not contiguous
    blocks=[]; start=None
    for i in range(len(mids)):
        hole = covered[i] and start is not None and i>0 and gap_after[idx[i-1]:idx[i]].any()
        brk = (not covered[i]) or (start is not None and i>0 and (mids[i]-mids[i-1])>max_phys_gap) or hole
        if brk and start is not None: blocks.append((start,i-1)); start=None
        if covered[i] and start is None: start=i
        if not covered[i]: start=None
    if start is not None: blocks.append((start,len(mids)-1))
    return blocks
