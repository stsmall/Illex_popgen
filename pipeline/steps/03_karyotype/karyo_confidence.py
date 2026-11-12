#!/usr/bin/env python
"""Posterior confidence of the chr2 karyotype calls (main-text Table 2, Supplementary Methods).

Refits the three-component Gaussian mixture of recluster_karyo.py (same input, settings and seed) to the local-PCA PC1
of all rows, checks that its labels reproduce karyotypes.tsv, and reports the share of the analysed individuals
(Baker et al. 2025 BAMs minus data/qc/excluded_14.tsv) whose largest posterior exceeds 0.95 and 0.99.
Writes data/qc/karyotype_posteriors.tsv.
"""
import os, subprocess
import numpy as np, pandas as pd
from sklearn.mixture import GaussianMixture

BASE = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype"
VCF = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/2/variants_filt.vcf.gz"
BCF = "/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"
BAKER = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/baker_633.txt"
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

df = pd.read_csv(f"{BASE}/chr2_karyo_coords.tsv", sep="\t")
df.insert(0, "sample", [s for s in subprocess.check_output([BCF, "query", "-l", VCF]).decode().split("\n") if s.strip()])
pc1 = df["PC1"].values.astype(float).reshape(-1, 1)
gm = GaussianMixture(n_components=3, n_init=10, random_state=0).fit(pc1)
order = np.argsort(gm.means_.ravel())
name = {order[0]: "AA", order[1]: "AB", order[2]: "BB"}
post = gm.predict_proba(pc1)
df["karyotype"] = [name[l] for l in post.argmax(1)]
df["max_posterior"] = post.max(1)

ref = pd.read_csv(f"{BASE}/karyotypes.tsv", sep="\t").set_index("sample")["karyotype"]
assert (df.set_index("sample")["karyotype"] == ref.reindex(df["sample"]).values).all(), "labels differ from karyotypes.tsv"

excl = set(pd.read_csv(f"{ROOT}/data/qc/excluded_14.tsv", sep="\t", comment="#")["sample"])
baker = set(open(BAKER).read().split())
d = df[df["sample"].isin(baker - excl)]
d[["sample", "karyotype", "max_posterior"]].to_csv(f"{ROOT}/data/qc/karyotype_posteriors.tsv", sep="\t", index=False)
print(f"n={len(d)} counts={d.karyotype.value_counts().to_dict()}")
for t in (0.9, 0.95, 0.99):
    print(f"max posterior > {t}: {(d.max_posterior > t).sum()} ({100 * (d.max_posterior > t).mean():.1f}%)")
print("not >0.95 by karyotype:", d[d.max_posterior <= 0.95].karyotype.value_counts().to_dict())
