#!/usr/bin/env python
# chr2 LD half-triangle from the saved r^2 matrix (no recompute). Robust method:
# rotate the upper-triangle matrix 45 deg with scipy.ndimage.rotate and crop to
# the top half. vmax tuned to the real signal (block mean ~0.02, flank ~0.006).
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import rotate

PRE = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/02_inv_pca/chr2_inv_ld"
r2 = np.load(PRE + "_r2.npy").astype(float)
pos = np.load(PRE + "_pos.npy").astype(float)
N = r2.shape[0]
r2 = np.clip(r2, 0, 1)                 # composite estimator can exceed 1; clip

# keep upper triangle (incl diagonal), NaN elsewhere so rotation background is transparent
M = np.full_like(r2, np.nan)
iu = np.triu_indices(N)
M[iu] = r2[iu]

# rotate 45 deg; nearest-neighbour to preserve values; NaN fill outside
Mrot = rotate(M, 45, reshape=True, order=0, cval=np.nan, prefilter=False)
h = Mrot.shape[0]
tri = Mrot[: h // 2 + 1, :]            # top half = the triangle sitting on the axis

cmap = plt.cm.magma.copy(); cmap.set_bad(alpha=0.0)   # NaN transparent
fig, ax = plt.subplots(figsize=(16, 5.5))
ax.set_facecolor("white")
im = ax.imshow(np.ma.masked_invalid(tri), aspect="auto", cmap=cmap,
               vmin=0, vmax=0.10, origin="upper", interpolation="nearest",
               extent=[pos.min()/1e6, pos.max()/1e6, 0, 1])
ax.set_yticks([])
ax.set_xlabel("chr2 position (Mb)")
ax.set_ylabel("SNP separation →")
ax.set_title("chr2 LD half-triangle (r²) — variants_filt, 2:55–85 Mb, 2500 SNPs (vmax=0.10)")
for bp in (60, 80):
    ax.axvline(bp, color="cyan", ls="--", lw=1.2, alpha=0.8)
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01); cb.set_label("r²")
fig.tight_layout()
OUT = PRE + "_triangle.png"
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print("wrote", OUT, "| block/flank mean r2 = 0.022/0.006")
