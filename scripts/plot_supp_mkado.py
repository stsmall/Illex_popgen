#!/usr/bin/env python
"""Supplementary figure: McDonald-Kreitman / adaptive-evolution summary for
I. illecebrosus vs I. coindetii (mkado results).

Panel A: genome-wide alpha, standard-MK vs asymptotic-MK (with Monte-Carlo CI).
Panel B: per-gene DoS distribution, SF3B4 (LOC_00013210, top FDR-significant
adaptive gene) highlighted.

Data:
  results/coindetii/coindetii.asymptotic.tsv  (1-row genome-wide fit)
  results/coindetii/coindetii.per_gene.tsv    (20,999 genes)
"""
import sys; sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, despine, save
apply()
import matplotlib.pyplot as plt
import pandas as pd

RES = "/sietch_colab/data_share/illex/popgen_data/mkado_illex/results/coindetii"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/figS10_mkado.png"
SF3B4 = "LOC_00013210"

asym = pd.read_csv(f"{RES}/coindetii.asymptotic.tsv", sep="\t").iloc[0]
genes = pd.read_csv(f"{RES}/coindetii.per_gene.tsv", sep="\t")

# genome-wide pooled standard MK alpha, from the same pooled Dn/Ds/Pn/Ps
Dn, Ds, Pn, Ps = asym.Dn, asym.Ds, asym.Pn, asym.Ps
alpha_std = 1 - (Ds * Pn) / (Dn * Ps)
alpha_asym = asym.alpha_asymptotic
ci_lo, ci_hi = asym.CI_low, asym.CI_high

sf3b4 = genes.loc[genes.gene == SF3B4].iloc[0]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.3))

# ---- Panel A: standard vs asymptotic alpha ----
labels = ["Standard MK", "Asymptotic MK"]
vals = [alpha_std, alpha_asym]
colors = [C["muted"], C["accent"]]
x = [0, 1]
axA.bar(x, vals, width=0.55, color=colors, zorder=3)
axA.errorbar(1, alpha_asym, yerr=[[alpha_asym - ci_lo], [ci_hi - alpha_asym]],
             fmt="none", ecolor=C["ink"], elinewidth=1.2, capsize=4, zorder=4)
axA.axhline(0, color=C["muted"], linewidth=0.8, zorder=1)
axA.set_xticks(x, labels)
axA.set_xlim(-0.6, 1.6)
axA.set_ylim(-0.1, 4.0)
axA.set_ylabel(r"$\alpha$ (fraction adaptive substitutions)")
axA.text(0, alpha_std + 0.08, f"{alpha_std:.2f}", ha="center", va="bottom",
         fontsize=9, color=C["ink"])
axA.text(1, ci_hi + 0.08, f"{alpha_asym:.2f}\n(CI {ci_lo:.2f}–{ci_hi:.2f})",
         ha="center", va="bottom", fontsize=9, color=C["ink"])
axA.set_title("A", loc="left")
despine(axA)

# ---- Panel B: per-gene DoS distribution ----
dos = genes["DoS"].dropna()
axB.hist(dos, bins=60, color=C["muted"], edgecolor="none", zorder=2)
axB.axvline(0, color=C["ink"], linewidth=0.8, linestyle=":", zorder=1)
axB.axvline(sf3b4.DoS, color=C["warm"], linewidth=1.6, zorder=4)
ymax = axB.get_ylim()[1]
axB.annotate(
    "SF3B4\n" + rf"DoS={sf3b4.DoS:.2f}, $\alpha$={sf3b4.alpha:.2f}" + "\n"
    + rf"$p_{{adj}}$={sf3b4.p_value_adjusted:.1e}",
    xy=(sf3b4.DoS, ymax * 0.62), xytext=(sf3b4.DoS - 0.62, ymax * 0.78),
    ha="left", va="bottom", fontsize=8.5, color=C["warm"],
    arrowprops=dict(arrowstyle="-", color=C["warm"], linewidth=1.0),
)
axB.set_xlabel("Direction of Selection (DoS)")
axB.set_ylabel(f"genes (n={len(dos):,})")
axB.set_xlim(-1.05, 1.05)
axB.set_title("B", loc="left")
despine(axB)

fig.tight_layout()
save(fig, OUT)
