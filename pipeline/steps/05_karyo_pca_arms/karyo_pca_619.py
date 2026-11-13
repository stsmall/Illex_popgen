#!/usr/bin/env python
"""Per-chromosome PCA of the thinned SNPs of the 619 analysed individuals (refit_pca_619.sh).

Same recipe as karyo_pca_panel.py: pg_gpu randomized PCA, 4 components, haplotype pairs averaged per
individual. Writes OUTDIR/karyo_pca_coords_619.tsv (chrom, sample, PC1, PC2), the input of
eta2_karyofix.py. Usage: karyo_pca_619.py OUTDIR
"""
import glob
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from pg_gpu import HaplotypeMatrix, decomposition

O = sys.argv[1]
BCF = "/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"


def key(f):
    c = os.path.basename(f).replace("chr", "").replace(".thin.vcf.gz", "")
    return (0, int(c)) if c.isdigit() else (1, c)


rows = []
for vcf in sorted(glob.glob(f"{O}/thinned/chr*.thin.vcf.gz"), key=key):
    c = os.path.basename(vcf).replace("chr", "").replace(".thin.vcf.gz", "")
    h = HaplotypeMatrix.from_vcf(vcf)
    samples = subprocess.check_output([BCF, "query", "-l", vcf]).decode().split()
    n = len(samples)
    coords, _ = decomposition.randomized_pca(h, n_components=4)
    coords = coords.get() if hasattr(coords, "get") else np.asarray(coords)
    if coords.shape[0] == 2 * n:
        coords = (coords[:n] + coords[n:2 * n]) / 2.0
    rows += [(c, s, a, b) for s, a, b in zip(samples, coords[:, 0], coords[:, 1])]
    print(f"chr{c}: {h.num_variants} SNPs, {n} individuals", flush=True)
pd.DataFrame(rows, columns=["chrom", "sample", "PC1", "PC2"]).to_csv(f"{O}/karyo_pca_coords_619.tsv", sep="\t", index=False)
print("wrote karyo_pca_coords_619.tsv")
