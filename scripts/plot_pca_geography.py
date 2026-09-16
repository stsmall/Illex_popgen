"""Supp: the autosomal-PCA spread does NOT track geography.
a. autosomal PCA (chr2 excluded) coloured by NAFO division.
b. same PCA coloured by latitude.
c. PC1 and PC2 vs latitude, with Pearson r (near zero) -> spread is not geographic.
"""
import os
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
_HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(_HERE, "..", "data")
from figstyle import apply, C, OKABE, despine
apply()
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from scipy.stats import pearsonr

SC = DATA
META = "/sietch_colab/data_share/illex/popgen_data/pggpu_illex/popstats/metadata.tsv"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/figS1_pca_geography.png"

ev = pd.read_csv(f"{SC}/auto_pca4.eigenvec", sep="\t")
ev = ev.rename(columns={"IID": "sample"})
meta = pd.read_csv(META, sep="\t")[["sample", "population", "latitude", "longitude"]]
d = ev.merge(meta, on="sample", how="inner")
print(f"joined {len(d)}/{len(ev)} PCA samples to metadata")

# clip extreme outliers for display (robust, matches Fig 1); keep for the correlation
for pc in ("PC1", "PC2"):
    med, mad = d[pc].median(), (d[pc] - d[pc].median()).abs().median()
    d = d[(d[pc] - med).abs() <= 8 * mad * 1.4826]

order = sorted(d["population"].unique(),
               key=lambda p: d[d["population"] == p]["latitude"].mean())

fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.0), gridspec_kw={"wspace": 0.42})

# a: colour by division
for i, p in enumerate(order):
    s = d[d["population"] == p]
    ax[0].scatter(s["PC1"], s["PC2"], s=12, color=OKABE[i % len(OKABE)], label=p,
                  alpha=0.8, linewidths=0)
ax[0].set_xlabel("PC1"); ax[0].set_ylabel("PC2")
ax[0].set_title("a   PCA by NAFO division", loc="left")
ax[0].legend(fontsize=7, ncol=2, handletextpad=0.2, columnspacing=0.8,
             borderaxespad=0.2, loc="best", markerscale=1.2)
despine(ax[0])

# b: colour by latitude
norm = mpl.colors.Normalize(d["latitude"].min(), d["latitude"].max())
sc = ax[1].scatter(d["PC1"], d["PC2"], s=12, c=d["latitude"], cmap="viridis",
                   norm=norm, alpha=0.85, linewidths=0)
ax[1].set_xlabel("PC1"); ax[1].set_ylabel("PC2")
ax[1].set_title("b   PCA by latitude", loc="left")
cb = fig.colorbar(sc, ax=ax[1], fraction=0.046, pad=0.03)
cb.set_label("latitude (°N)", fontsize=8); cb.ax.tick_params(labelsize=7)
despine(ax[1])

# c: PC vs latitude, with r
r1, p1 = pearsonr(d["latitude"], d["PC1"])
r2, p2 = pearsonr(d["latitude"], d["PC2"])
ax[2].scatter(d["latitude"], d["PC1"], s=10, color=C["accent"], alpha=0.6,
              linewidths=0, label=f"PC1  (r={r1:+.2f})")
ax[2].scatter(d["latitude"], d["PC2"], s=10, color=C["warm"], alpha=0.6,
              linewidths=0, label=f"PC2  (r={r2:+.2f})")
ax[2].set_xlabel("latitude (°N)"); ax[2].set_ylabel("PC score")
ax[2].set_title("c   PC vs latitude", loc="left")
ax[2].legend(fontsize=8, loc="best")
despine(ax[2])

fig.savefig(OUT, dpi=220, bbox_inches="tight")
print("wrote", OUT)
print(f"PC1~lat r={r1:+.3f} (p={p1:.2f}); PC2~lat r={r2:+.3f} (p={p2:.2f})")
