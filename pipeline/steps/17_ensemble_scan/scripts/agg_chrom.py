#!/usr/bin/env python
"""Aggregate one chromosome's 45 plink2 .fst.var files to per-SNP mean FST.
Same filters/clip as the simulation null:
  - per-pair: OBS_CT>=40, drop nan
  - clip per-pair WC_FST to [0,1]
  - per-SNP mean across pairs, require >=20 informative pairs
Writes <outdir>/perchrom_<chrom>.tsv  (CHROM POS meanf maxf npair)
"""
import sys, glob
import numpy as np, pandas as pd

chrom = sys.argv[1]           # e.g. "1" or "Z"
d     = sys.argv[2]           # dir containing fst_<chrom>.*.fst.var
prefix = f"{d}/fst_{chrom}."
files = sorted(glob.glob(prefix + "*.fst.var"))
if not files:
    sys.exit(f"no fst.var for chrom {chrom} in {d}")
frames = []
for f in files:
    x = pd.read_csv(f, sep="\t", usecols=[1,3,4], names=["POS","obs","fst"], header=0)
    x = x[x["obs"] >= 40].dropna(subset=["fst"])
    if len(x) == 0:
        continue
    x["fst"] = np.clip(x["fst"].to_numpy(), 0.0, 1.0)
    frames.append(x[["POS","fst"]])
allf = pd.concat(frames, ignore_index=True)
g = allf.groupby("POS", sort=True).agg(
        meanf=("fst", "mean"), maxf=("fst", "max"), npair=("fst", "size")).reset_index()
g = g[g["npair"] >= 20].reset_index(drop=True)
g.insert(0, "CHROM", chrom)
g.to_csv(f"{d}/perchrom_{chrom}.tsv", sep="\t", index=False)
print(f"chrom {chrom}: {len(g)} per-SNP records (n_pairfiles={len(files)})", flush=True)
