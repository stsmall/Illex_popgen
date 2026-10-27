#!/usr/bin/env python
"""Per-gene NG86 synonymous and nonsynonymous divergence (Illex vs Sthenoteuthis) from MACSE codon
alignments, then the Z/autosome ratio R and the male:female mutation ratio alpha = (3R-2)/(4-3R)
(Miyata et al. 1987), as in Coffing et al. 2025. Codons with a gap, frameshift (!), N or a stop in
either sequence are dropped. Genes in the 22 anomalous-diversity blocks and the chr2 inversion are
excluded. R is computed from mean per-gene dS (as in Coffing et al.) and from pooled synonymous
differences over pooled synonymous sites (JC-corrected). CIs: 10,000 bootstrap resamples of genes.
Usage: ds_ng86.py PAIRS_DIR BLOCKS_BED OUT_PREFIX"""
import sys, glob, os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from itertools import permutations, product
P, BLK, OUT = sys.argv[1:4]
# standard genetic code, Nei & Gojobori (1986) site and difference counting
B = "TCAG"; AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
code = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
stops = {c for c, a in code.items() if a == "*"}
def syn_sites(c):
    s = 0.0
    for i in range(3):
        for b in "ACGT":
            if b == c[i]: continue
            m = c[:i] + b + c[i+1:]
            if code[m] != "*" and code[m] == code[c]: s += 1 / 3
    return s
SS = {c: syn_sites(c) for c in code if c not in stops}
def codon_diff(u, v):
    pos = [i for i in range(3) if u[i] != v[i]]
    if not pos: return 0.0, 0.0
    sd = nd = 0.0; npath = 0
    for order in permutations(pos):
        cur, s_, n_, ok = u, 0, 0, True
        for i in order:
            nxt = cur[:i] + v[i] + cur[i+1:]
            if code[nxt] == "*": ok = False; break
            if code[nxt] == code[cur]: s_ += 1
            else: n_ += 1
            cur = nxt
        if ok: sd += s_; nd += n_; npath += 1
    return (sd / npath, nd / npath) if npath else (np.nan, np.nan)
def jc(p): return -0.75 * np.log(1 - 4 * p / 3) if p < 0.75 else np.nan
def readfa(f):
    d, k = {}, None
    for l in open(f):
        if l.startswith(">"): k = l[1:].strip(); d[k] = []
        else: d[k].append(l.strip())
    return {k: "".join(v).upper() for k, v in d.items()}
pairs = pd.read_csv(f"{P}/pairs.tsv", sep="\t", dtype={"chrom": str})
blk = pd.read_csv(BLK, sep="\t", header=None, names=["c", "s", "e"], dtype={"c": str})
rows = []
for r in pairs.itertuples():
    f = f"{P}/{r.tx}.NT.fa"
    if not os.path.exists(f): continue
    a = readfa(f)
    if "illex" not in a or "stheno" not in a: continue
    x, y = a["illex"].replace("!", "-"), a["stheno"].replace("!", "-")
    c1, c2 = [], []
    for i in range(0, min(len(x), len(y)) - 2, 3):
        u, v = x[i:i+3], y[i:i+3]
        if set(u + v) - set("ACGT") or u in stops or v in stops: continue
        c1.append(u); c2.append(v)
    if len(c1) < 50: continue
    S = sum((SS[u] + SS[v]) / 2 for u, v in zip(c1, c2)); N = 3 * len(c1) - S
    dd = np.array([codon_diff(u, v) for u, v in zip(c1, c2)]); Sd, Nd = np.nansum(dd[:, 0]), np.nansum(dd[:, 1])
    rows.append({"tx": r.tx, "chrom": r.chrom, "start": r.start, "end": r.end, "codons": len(c1),
                 "S": S, "N": N, "Sd": Sd, "Nd": Nd, "dS": jc(Sd / S), "dN": jc(Nd / N)})
g = pd.DataFrame(rows)
inv = (g.chrom == "2") & (g.end >= 60_540_000) & (g.start <= 79_500_000)
inb = np.zeros(len(g), bool)
for b in blk.itertuples():
    inb |= (g.chrom == b.c) & (g.end >= b.s) & (g.start <= b.e)
g["excluded"] = inv | inb
g.to_csv(f"{OUT}.per_gene.tsv", sep="\t", index=False)
k = g[~g.excluded & g.dS.notna()]
Z, A = k[k.chrom == "Z"], k[k.chrom != "Z"]
alpha = lambda R: (3 * R - 2) / (4 - 3 * R)
def stats(z, a):
    Rm = z.dS.mean() / a.dS.mean()
    Rp = jc(z.Sd.sum() / z.S.sum()) / jc(a.Sd.sum() / a.S.sum())
    return Rm, Rp
Rm, Rp = stats(Z, A)
rng = np.random.default_rng(1); bm, bp = [], []
for _ in range(10000):
    z = Z.iloc[rng.integers(0, len(Z), len(Z))]; a = A.iloc[rng.integers(0, len(A), len(A))]
    m, p = stats(z, a); bm.append(m); bp.append(p)
ci = lambda v: np.percentile(v, [2.5, 97.5])
res = {"n_Z": len(Z), "n_A": len(A), "mean_dS_Z": Z.dS.mean(), "mean_dS_A": A.dS.mean(),
       "median_dS_Z": Z.dS.median(), "median_dS_A": A.dS.median(),
       "mean_dN_Z": Z.dN.mean(), "mean_dN_A": A.dN.mean(),
       "R_mean": Rm, "R_mean_CI": ci(bm).tolist(), "alpha_mean": alpha(Rm), "alpha_mean_CI": alpha(ci(bm)).tolist(),
       "R_pooled": Rp, "R_pooled_CI": ci(bp).tolist(), "alpha_pooled": alpha(Rp), "alpha_pooled_CI": alpha(ci(bp)).tolist()}
from scipy.stats import mannwhitneyu
res["MW_dS_p_two_sided"] = float(mannwhitneyu(Z.dS, A.dS).pvalue)
res["MW_dN_p_two_sided"] = float(mannwhitneyu(Z.dN.dropna(), A.dN.dropna()).pvalue)
import json; json.dump(res, open(f"{OUT}.summary.json", "w"), indent=1, default=float)
for kk, v in res.items(): print(kk, v)
