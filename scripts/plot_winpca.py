"""Supp Fig S1: windowed PC1 along chr2 around the inversion, by karyotype.

Inputs (analysis/steps/05_karyo_pca_arms/refit_619, from pipeline/steps/05_karyo_pca_arms/refit_pca_619.sh):
chr2_winpca_619_coords.npy (windows x 619 samples x 10 PCs, sample order chr2_619_samples.txt),
chr2_winpca_619_windows.tsv, and data/qc/karyotypes_619.tsv (karyotype_final). PC1 sign is arbitrary per window, so each window is
oriented so that 2+/2+ individuals have positive mean PC1. 18 random individuals per karyotype.
"""
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

R = Path("/sietch_colab/data_share/illex/popgen_data/analysis/steps/05_karyo_pca_arms/refit_619")
OUT = Path(__file__).resolve().parent.parent / "supp_figures" / "figS1_winpca.png"
BP_L, BP_R = 60.54, 79.50

coords = np.load(R / "chr2_winpca_619_coords.npy")[:, :, 0]   # windows x samples
win = pd.read_csv(R / "chr2_winpca_619_windows.tsv", sep="\t")
names = (R / "chr2_619_samples.txt").read_text().split()
assert len(names) == coords.shape[1] == 619
kar = pd.read_csv(Path(__file__).resolve().parents[1] / "data/qc/karyotypes_619.tsv", sep="\t", comment="#").set_index("sample")["karyotype_final"]
lab = np.array([kar[n] for n in names])

sign = np.sign(np.nanmean(coords[:, lab == "BB"], axis=1))
sign[sign == 0] = 1
pc1 = coords * sign[:, None]
x = win["center"].to_numpy() / 1e6

rng = np.random.default_rng(1)
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.axvspan(BP_L, BP_R, color="#d9d9d9", lw=0, zorder=0)
for key, disp, col in (("AA", "2a/2a", "#2f78d0"), ("AB", "2a/2+", "#e8643a"), ("BB", "2+/2+", "#1fa77a")):
    idx = rng.choice(np.flatnonzero(lab == key), 18, replace=False)
    for k, i in enumerate(idx):
        ax.plot(x, pc1[:, i], lw=0.9, color=col, alpha=0.5, label=f"{disp} (n=18 shown)" if k == 0 else None)
ax.set_xlim(58, 82)
ax.set_xlabel("chr2 position (Mb)")
ax.set_ylabel("PC1 (window, sign-oriented)")
ax.legend(frameon=False, loc="upper right")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print("wrote", OUT, {k: int((lab == k).sum()) for k in ("AA", "AB", "BB")})
