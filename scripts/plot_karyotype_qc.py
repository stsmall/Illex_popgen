#!/usr/bin/env python
"""Supplementary figure: sequencing depth, call rate, called heterozygosity and sampling by chr2 karyotype group.

Inputs: analysis/steps/03_karyotype/m4_qc/perind_chr20_maf05.tsv (written by
pipeline/steps/03_karyotype/m4_qc/karyotype_qc.py: per-individual autosomal read count from samtools idxstats,
chr20 call rate and Ho/He at biallelic PASS SNPs with MAF >= 0.05 from plink2, final karyotype and division).
Panels: (a) autosomal reads, (b) chr20 call rate, (c) Ho/He against reads, (d) karyotype composition per division.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import apply, C, KARYO, despine, save  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

apply()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = pd.read_csv("/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/m4_qc/perind_chr20_maf05.tsv",
                sep="\t", index_col=0)
DISP = {"AA": "2a/2a", "AB": "2a/2+", "BB": "2+/2+"}
ORDER = ["AA", "AB", "BB"]
rng = np.random.default_rng(1)

fig, ax = plt.subplots(2, 2, figsize=(7.2, 6.0))
for a, col, lab in ((ax[0, 0], "reads", "autosomal reads (millions)"), (ax[0, 1], "call", "chr20 call rate")):
    for i, k in enumerate(ORDER):
        v = D.loc[D.K == k, col].dropna().to_numpy()
        a.scatter(i + rng.uniform(-0.18, 0.18, v.size), v, s=6, color=KARYO[k], alpha=0.55, lw=0)
        a.plot([i - 0.28, i + 0.28], [np.median(v)] * 2, color=C["ink"], lw=1.6)
    a.set_xticks(range(3), [f"{DISP[k]}\n(n={int((D.K == k).sum())})" for k in ORDER])
    a.set_ylabel(lab)
    despine(a)
a = ax[1, 0]
for k in ORDER:
    s = D[D.K == k]
    a.scatter(s.reads, s.hohe, s=6, color=KARYO[k], alpha=0.6, lw=0, label=DISP[k])
a.set_xscale("log")
a.set_xlabel("autosomal reads (millions)")
a.set_ylabel("called Ho/He, chr20")
a.legend(frameon=False, fontsize=8, markerscale=2)
despine(a)
a = ax[1, 1]
ct = pd.crosstab(D["div"], D.K)
ct = ct[ct.sum(1) >= 5]
frac = ct.div(ct.sum(1), axis=0)
left = np.zeros(len(frac))
for k in ORDER:
    a.barh(range(len(frac)), frac[k], left=left, color=KARYO[k], height=0.75, edgecolor="white", lw=0.6, label=DISP[k])
    left += frac[k].to_numpy()
a.set_yticks(range(len(frac)), [f"{d} ({int(n)})" for d, n in zip(frac.index, ct.sum(1))])
a.set_xlabel("karyotype share")
a.set_xlim(0, 1)
despine(a)
for a, l in zip(ax.ravel(), "abcd"):
    a.set_title(l, loc="left", fontweight="bold")
fig.tight_layout()
save(fig, os.path.join(ROOT, "supp_figures", "figS27_karyotype_qc.png"))
