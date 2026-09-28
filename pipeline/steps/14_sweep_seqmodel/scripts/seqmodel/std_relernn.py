import numpy as np, pandas as pd
_COL = {"corrected":"recombRate","lower":"CI95LO","upper":"CI95HI"}

def standardize(path, which="corrected"):
    df = pd.read_csv(path, sep="\t")
    df["chrom"] = df["chrom"].astype(str).str.replace(r"^b'|'$","",regex=True)
    if ((df.CI95LO > df.recombRate) | (df.recombRate > df.CI95HI)).any():
        raise ValueError("CI ordering violated: need CI95LO <= recombRate <= CI95HI")
    out = df[["chrom","start","end"]].copy()
    out["rate"] = df[_COL[which]].astype(float).clip(lower=0.0)
    out = out.sort_values(["chrom","start"]).reset_index(drop=True)
    for c,g in out.groupby("chrom"):
        s,e = g.start.to_numpy(), g.end.to_numpy()
        if (e <= s).any(): raise ValueError(f"{c}: end<=start")
        if (s[1:] < e[:-1]).any(): raise ValueError(f"{c}: overlapping intervals")
    return out

def integrated_length_morgans(df, chrom):
    g = df[df.chrom == str(chrom)]
    return float(((g.end - g.start) * g.rate).sum())
