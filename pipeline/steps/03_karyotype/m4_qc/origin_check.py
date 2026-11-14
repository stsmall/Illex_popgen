#!/usr/bin/env python
"""Single or multiple origin of 2a (review M2): genotype distances among 2a/2a homozygotes near the breakpoints.

If 2a was founded more than once, 2a/2a individuals in the first and last megabase of the inversion should fall into
deeply divergent groups. For each edge (chr2:60.54-61.54 and 78.50-79.50 Mb), biallelic PASS SNPs with minor-allele
frequency >= 0.05 among the 2a/2a and 2+/2+ individuals are read from the callset; the distance between two
individuals is the mean absolute dosage difference / 2 over sites called in both. Reports within-2a/2a, within-2+/2+
and between-group distances, the leading principal component of the 2a/2a individuals alone, and a Gaussian-mixture
test (BIC, k = 1 vs 2) of that component. SQ233K026 is left out (excluded sample).
"""
import subprocess

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

BCF = "/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"
VCF = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/2/variants_filt.vcf.gz"
K = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/illex_manuscript/data/qc/karyotypes_619.tsv", sep="\t",
                comment="#").set_index("sample").karyotype_final.drop("SQ233K026")
AA, BB = K[K == "AA"].index.tolist(), K[K == "BB"].index.tolist()


def dosages(region, samples):
    q = subprocess.run([BCF, "query", "-r", region, "-s", ",".join(samples), "-f", "[%GT\t]\n", "-i",
                        'FILTER="PASS" && TYPE="snp" && N_ALT=1', VCF], capture_output=True, text=True, check=True)
    code = {"0/0": 0, "0|0": 0, "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}
    g = np.array([[code.get(x, np.nan) for x in l.rstrip("\t").split("\t")] for l in q.stdout.splitlines()], float)
    return g


def dist(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return np.abs(a[ok] - b[ok]).mean() / 2 if ok.sum() > 100 else np.nan


for name, region in (("left edge", "2:60540000-61540000"), ("right edge", "2:78500000-79500000")):
    g = dosages(region, AA + BB)
    af = np.nanmean(g, 1) / 2
    g = g[(np.minimum(af, 1 - af) >= 0.05)]
    ga, gb = g[:, :len(AA)], g[:, len(AA):]
    rng = np.random.default_rng(1)
    pa = [dist(ga[:, i], ga[:, j]) for i, j in rng.integers(0, len(AA), (3000, 2)) if i != j]
    pb = [dist(gb[:, i], gb[:, j]) for i, j in rng.integers(0, len(BB), (3000, 2)) if i != j]
    pab = [dist(ga[:, i], gb[:, j]) for i, j in zip(rng.integers(0, len(AA), 3000), rng.integers(0, len(BB), 3000))]
    x = ga.copy()
    col = np.nanmean(x, 1, keepdims=True)
    x = np.where(np.isfinite(x), x, col) - col
    u, sv, vt = np.linalg.svd(x, full_matrices=False)
    pc1 = vt[0] * sv[0]
    X = pc1.reshape(-1, 1)
    bic = {k: GaussianMixture(k, n_init=4, random_state=0).fit(X).bic(X) for k in (1, 2)}
    print(f"{name}: {g.shape[0]} SNPs; distance within 2a/2a median {np.nanmedian(pa):.4f} (99th pct "
          f"{np.nanpercentile(pa, 99):.4f}), within 2+/2+ {np.nanmedian(pb):.4f}, between {np.nanmedian(pab):.4f}; "
          f"2a/2a PC1 variance share {sv[0] ** 2 / (sv ** 2).sum():.3f}, dBIC(k1-k2) {bic[1] - bic[2]:.1f}")
