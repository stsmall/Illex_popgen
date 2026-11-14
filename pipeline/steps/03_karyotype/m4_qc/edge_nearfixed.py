#!/usr/bin/env python
"""Near-fixed differences by distance from the nearer breakpoint (review minor 14).

SNPs variable among the 223 2a/2a and 97 2+/2+ individuals between the breakpoints (60.54-79.50 Mb), from the karyofix
allele-frequency tables (positions duplicated by multiallelic records dropped). Near-fixed = |AF_2a - AF_2+| >= 0.9.
Reports, per 1-Mb bin of distance from the nearer breakpoint, the SNP count, near-fixed count and near-fixed share per
SNP, the per-Mb densities of the outermost megabase against the interior, and the two breakpoints separately.
Writes analysis/steps/03_karyotype/m4_qc/edge/nearfixed_by_distance.tsv.
"""
import numpy as np
import pandas as pd

K = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/karyofix/"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/m4_qc/edge/nearfixed_by_distance.tsv"
L, R = 60_540_000, 79_500_000


def af(f):
    d = pd.read_csv(K + f, sep="\t", header=None, names=["c", "pos", "af"])
    d = d[(d.pos >= L) & (d.pos <= R)]
    return d.drop_duplicates("pos", keep=False).set_index("pos").af


d = pd.concat([af("AA223_chr2_10_80_af.txt").rename("fA"), af("BB97_chr2_10_80_af.txt").rename("fB")],
              axis=1, join="inner").dropna()
d = d[((d.fA > 0) & (d.fA < 1)) | ((d.fB > 0) & (d.fB < 1)) | (d.fA != d.fB)]
d["dist"] = np.minimum(d.index - L, R - d.index) / 1e6
d["nf"] = (d.fA - d.fB).abs() >= 0.9
d["bin"] = np.floor(d.dist).clip(upper=9).astype(int)
g = d.groupby("bin").agg(snps=("nf", "size"), nearfixed=("nf", "sum"))
g["share"] = g.nearfixed / g.snps
g.to_csv(OUT, sep="\t")
print(g.to_string())
outer, inner = g.loc[0], g.loc[1:].sum()
inner_mb = (R - L) / 1e6 - 2.0
print(f"outer 1 Mb (2 Mb total): share {outer.nearfixed / outer.snps:.4f}, {outer.nearfixed / 2:.0f} near-fixed and "
      f"{outer.snps / 2:.0f} SNPs per Mb")
print(f"interior ({inner_mb:.2f} Mb): share {inner.nearfixed / inner.snps:.4f}, {inner.nearfixed / inner_mb:.0f} near-fixed "
      f"and {inner.snps / inner_mb:.0f} SNPs per Mb")
for side, m in (("left", d.index - L < 1e6), ("right", R - d.index < 1e6)):
    print(f"{side} outer Mb: {int(d[m].nf.sum())}/{int(m.sum())} = {d[m].nf.mean():.4f}")
