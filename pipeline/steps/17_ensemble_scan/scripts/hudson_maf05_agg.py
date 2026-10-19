#!/usr/bin/env python
"""Per-SNP geographic Hudson FST across the 45 NAFO-division pairs (one chromosome).

Input: plink2 --freq counts --loop-cats POP output for one chromosome, i.e. one
<prefix>.<DIV>.acount per division (ALT_CTS, OBS_CT = allele observations), on SNPs
already filtered to MAF > 0.05 in the pooled 632 individuals.

Per pair (Bhatia et al. 2013, Genome Res 23:1514), with allele frequencies p1, p2 and
allele sample sizes n1, n2:
    num = (p1 - p2)^2 - p1(1 - p1)/(n1 - 1) - p2(1 - p2)/(n2 - 1)
    den = p1(1 - p2) + p2(1 - p1)
A pair is informative at a SNP if it has >= 40 genotyped individuals in total
(n1/2 + n2/2 >= 40), at least 2 alleles in each division, and den > 0.
Per SNP: ratio of averages over informative pairs, sum(num) / sum(den), kept if
>= 20 of the 45 pairs are informative. Values are not clipped.
Writes <outdir>/hudson_<chrom>.tsv: CHROM POS fst npair p_pooled.
With a 4th argument "npy", writes <outdir>/<chrom>.fst.npy instead: float32 FST for every
MAF-filtered SNP in input order, NaN where fewer than 20 pairs are informative (used for the
permutation null, so all permutations align SNP by SNP), plus <outdir>/<chrom>.pos.npy.
"""
import sys, glob, itertools
import numpy as np, pandas as pd

chrom, prefix, out = sys.argv[1], sys.argv[2], sys.argv[3]
mode = sys.argv[4] if len(sys.argv) > 4 else "tsv"
files = sorted(glob.glob(f"{prefix}.D*.acount"))
assert len(files) == 10, f"expected 10 divisions, got {len(files)}"
alt, obs, pos = [], [], None
for f in files:
    x = pd.read_csv(f, sep="\t", usecols=["ID", "ALT_CTS", "OBS_CT"])
    if pos is None:
        ids = x["ID"].to_numpy()
        pos = np.array([int(s.split(":")[1]) for s in ids])
    else:
        assert (x["ID"].to_numpy() == ids).all()
    alt.append(x["ALT_CTS"].to_numpy(np.float64)); obs.append(x["OBS_CT"].to_numpy(np.float64))
A = np.vstack(alt); O = np.vstack(obs)
with np.errstate(divide="ignore", invalid="ignore"):
    P = np.where(O > 0, A / O, np.nan)
    num_sum = np.zeros(A.shape[1]); den_sum = np.zeros(A.shape[1]); npair = np.zeros(A.shape[1], int)
    for i, j in itertools.combinations(range(10), 2):
        p1, p2, n1, n2 = P[i], P[j], O[i], O[j]
        num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        den = p1 * (1 - p2) + p2 * (1 - p1)
        ok = (n1 / 2 + n2 / 2 >= 40) & (n1 >= 2) & (n2 >= 2) & np.isfinite(num) & (den > 0)
        num_sum += np.where(ok, num, 0.0); den_sum += np.where(ok, den, 0.0); npair += ok
    keep = npair >= 20
    fst = num_sum[keep] / den_sum[keep]
if mode == "npy":
    full = np.full(A.shape[1], np.nan, np.float32); full[keep] = fst
    np.save(f"{out}/{chrom}.fst.npy", full)
    if not __import__("os").path.exists(f"{out}/../{chrom}.pos.npy"):
        np.save(f"{out}/../{chrom}.pos.npy", pos.astype(np.int64))
        with np.errstate(divide="ignore", invalid="ignore"):
            np.save(f"{out}/../{chrom}.maf.npy", np.minimum(A.sum(0) / O.sum(0), 1 - A.sum(0) / O.sum(0)).astype(np.float32))
    sys.exit(0)
with np.errstate(divide="ignore", invalid="ignore"):
    p_pooled = A.sum(0) / O.sum(0)
pd.DataFrame({"CHROM": chrom, "POS": pos[keep], "fst": fst.round(6), "npair": npair[keep],
              "p_pooled": p_pooled[keep].round(4)}).to_csv(f"{out}/hudson_{chrom}.tsv", sep="\t", index=False)
print(f"chrom {chrom}: {A.shape[1]} MAF>0.05 SNPs, {int(keep.sum())} with >=20 informative pairs", flush=True)
