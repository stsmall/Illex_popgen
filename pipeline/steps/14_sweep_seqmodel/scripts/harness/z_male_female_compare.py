#!/usr/bin/env python3
"""Compare chrZ diversity (pi, thetaW, Tajima's D) between the 330 males (Z-diploid, the diploSHIC
scan sample) and the 289 females (ZW = Z-HAPLOID). A real Z sweep is shared by both sexes (same Z
pool) -> female pattern should track male; female Z heterozygosity should be ~0 (hemizygous)."""
import sys
import numpy as np
import allel

MALE = sys.argv[1]; FEM = sys.argv[2]; WIN = 100000; ZLEN = 38234284

def load(path):
    c = allel.read_vcf(path, fields=['variants/POS', 'calldata/GT'])
    pos = c['variants/POS']; gt = allel.GenotypeArray(c['calldata/GT'])
    return pos, gt

def winstats(pos, gt):
    ac = gt.count_alleles()
    pi, win, nb, _ = allel.windowed_diversity(pos, ac, size=WIN, start=1, stop=ZLEN)
    tw, _, _, _   = allel.windowed_watterson_theta(pos, ac, size=WIN, start=1, stop=ZLEN)
    d, _, _       = allel.windowed_tajima_d(pos, ac, size=WIN, start=1, stop=ZLEN)
    # mean per-genotype heterozygosity per window (hemizygosity check)
    het_site = gt.count_het(axis=1) / gt.n_samples
    hetw = np.full(len(win), np.nan)
    for i,(a,b) in enumerate(win):
        m = (pos>=a)&(pos<=b)
        if m.sum(): hetw[i] = het_site[m].mean()
    return win, pi, tw, d, hetw

print("loading male...", file=sys.stderr); mp, mg = load(MALE)
print("loading female...", file=sys.stderr); fp, fg = load(FEM)
print(f"males n={mg.n_samples} sites={mg.n_variants} | females n={fg.n_samples} sites={fg.n_variants}")
mw, mpi, mtw, md, mhet = winstats(mp, mg)
fw, fpi, ftw, fd, fhet = winstats(fp, fg)

# align windows (same grid)
ok = np.isfinite(mpi) & np.isfinite(fpi) & np.isfinite(md) & np.isfinite(fd)
def corr(x,y):
    xy=np.isfinite(x)&np.isfinite(y); return np.corrcoef(x[xy],y[xy])[0,1]
print("\n=== genome-wide chrZ, male vs female (per-100kb window) ===")
print(f"windows compared: {ok.sum()}")
print(f"  pi:      male median {np.nanmedian(mpi):.5f}  female median {np.nanmedian(fpi):.5f}  corr r={corr(mpi,fpi):.3f}")
print(f"  thetaW:  male median {np.nanmedian(mtw):.5f}  female median {np.nanmedian(ftw):.5f}  corr r={corr(mtw,ftw):.3f}")
print(f"  TajD:    male median {np.nanmedian(md):.3f}   female median {np.nanmedian(fd):.3f}   corr r={corr(md,fd):.3f}")
print(f"  female Z heterozygosity (mean per-window): median {np.nanmedian(fhet):.4f}  (near 0 => hemizygous, Z-linked OK)")
print(f"  male   Z heterozygosity: median {np.nanmedian(mhet):.4f}")

# candidate windows
print("\n=== chrZ focal candidates: male vs female ===")
cands = [(12700001,"hard0.95"),(20800001,"hard0.58"),(20900001,"hard0.58"),(21600001,"hard0.95"),
         (23600001,"hard0.97"),(24500001,"soft0.92"),(25400001,"soft0.46"),(30700001,"soft0.78"),(31400001,"hard0.63")]
wid = {int(a):i for i,(a,b) in enumerate(mw)}
print(f"{'window':<14}{'class':<10}{'pi_M':>8}{'pi_F':>8}{'tW_M':>8}{'tW_F':>8}{'D_M':>7}{'D_F':>7}{'hetF':>7}")
for start,lab in cands:
    i = wid.get(start)
    if i is None: continue
    print(f"Z:{start/1e6:<11.2f}{lab:<10}{mpi[i]:>8.5f}{fpi[i]:>8.5f}{mtw[i]:>8.5f}{ftw[i]:>8.5f}{md[i]:>7.2f}{fd[i]:>7.2f}{fhet[i]:>7.3f}")
