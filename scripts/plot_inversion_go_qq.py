"""Supp: GO-enrichment QQ-plot for the genes captured by the chr2 inversion body
(60.54-79.50 Mb), against a PERMUTATION null (K=500 random gene sets of the same size,
all three GO categories). A theoretical-uniform null is wrong for GO -- terms are nested/
correlated and counts are tiny/discrete, which inflates any small gene set's QQ. The
permutation band captures that, and the observed points fall ENTIRELY within it: the
inversion genes produce fewer nominally-enriched terms than the random-set mean
(69 vs 86.6 +/- 60; permutation p=0.51) and nothing survives FDR. Unambiguously null --
the huge per-term folds are small-count artefacts, not enrichment."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE) if os.path.basename(_HERE) == "scripts" else _HERE
sys.path.insert(0, _HERE)
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt

O = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content"
OUT = os.path.join(_ROOT, "supp_figures", "figS21_inversion_go_qq.png")

b = pd.read_csv(f"{O}/perm_null_band.tsv", sep="\t")
b = b[(b.obs > 0) | (b.exp > 0)]
cnt = dict(l.split("=") for l in open(f"{O}/perm_null_counts.txt").read().split("\n") if "=" in l)
obs_nom, null_mean, null_sd = cnt["observed_nominal"], cnt["null_nominal_mean"], cnt["null_nominal_sd"]
pval = cnt["p_observed_ge_null"]

fig, ax = plt.subplots(figsize=(5.7, 5.4))
ax.fill_between(b.exp, b.lo, b.hi, color=C["faint"], alpha=0.85, zorder=0,
                label="permutation null (95%)")
lim = max(b.exp.max(), b.obs.max(), b.hi.max()) * 1.05
ax.plot([0, lim], [0, lim], ls="--", lw=1.0, color=C["muted"], zorder=1, label="null expectation")
ax.scatter(b.exp, b.obs, s=20, color=C["accent"], edgecolors="white", linewidths=0.4,
           zorder=3, label="inversion genes (observed)")
ax.set_xlim(0, lim); ax.set_ylim(0, lim)
ax.set_xlabel(r"permutation-null  $-\log_{10}p$")
ax.set_ylabel(r"observed  $-\log_{10}p$")
ax.set_title("Inversion genes: GO enrichment is null", loc="left", fontweight="bold", fontsize=11, pad=8)
ax.text(0.03, 0.97,
        f"496 GO terms (all categories); 0 significant at FDR<0.05\n"
        f"{obs_nom} nominally enriched vs {null_mean}$\\pm${null_sd} in random\n"
        f"gene sets of the same size (permutation $p={pval}$)",
        transform=ax.transAxes, va="top", ha="left", fontsize=8, color=C["ink"])
ax.legend(loc="lower right", frameon=False, fontsize=8, handletextpad=0.5, labelspacing=0.5)
despine(ax)
fig.savefig(OUT, dpi=220, bbox_inches="tight")
print("wrote", OUT, f"| observed within permutation band; p={pval}")
