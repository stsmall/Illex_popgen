"""LD triangles (genotype r^2, common SNPs, binned) for the chr2 inversion (positive control) and the
strongest elevated-diversity blocks; each spans the region +/- 1/3 of its length of flank. Colour scale shared."""
import os, sys
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
_HERE = os.path.dirname(os.path.abspath(__file__)); _ROOT = os.path.dirname(_HERE); sys.path.insert(0, _HERE)
from figstyle import apply, C
apply()
import numpy as np, matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_perpop/ld"
OUT = os.path.join(_ROOT, "supp_figures", "figS24_block_ld_triangles.png")
REG = [("chr2inv_control", "chr2 inversion (positive control)"), ("chr1_24_31", "chr1 block"),
       ("chr30_17_18", "chr30 block (6.1× π)"), ("chr23_0_6", "chr23 block (4.3× π)"),
       ("chr38_61_68", "chr38 block (4.1× π)"), ("chr35_60_65", "chr35 block (3.4× π)")]
NB = 120; VMAX = 0.008
fig, axs = plt.subplots(2, 3, figsize=(13, 6.4)); axs = axs.ravel()
for ax, (n, title) in zip(axs, REG):
    z = np.load(f"{D}/{n}.npz"); pos, r2, bs, be = z["pos"], z["r2"].astype(float), int(z["bs"]), int(z["be"])
    lo, hi = pos.min(), pos.max(); edges = np.linspace(lo, hi + 1, NB + 1); b = np.digitize(pos, edges) - 1
    M = np.full((NB, NB), np.nan)
    np.fill_diagonal(r2, np.nan)
    for i in range(NB):
        ii = b == i
        if not ii.any(): continue
        for j in range(i, NB):
            jj = b == j
            if jj.any(): M[i, j] = np.nanmean(r2[np.ix_(ii, jj)]) if j > i else np.nan
    e = np.arange(NB + 1); I, J = np.meshgrid(e, e, indexing="ij")
    X = (I + J) / 2.0; Y = (J - I) / 2.0
    im = ax.pcolormesh(X, Y, M, cmap="magma_r", vmin=0, vmax=VMAX, shading="flat", rasterized=True)
    ax.set_xlim(0, NB); ax.set_ylim(-4, NB / 2 + 1); ax.set_aspect("equal"); ax.axis("off")
    x0, x1 = (bs - lo) / (hi - lo) * NB, (be - lo) / (hi - lo) * NB
    ax.plot([x0, x1], [-2, -2], color=C["warm"], lw=4, solid_capstyle="butt")
    ax.text(0, -3.5, f"{lo/1e6:.1f}", fontsize=7, ha="left", va="top"); ax.text(NB, -3.5, f"{hi/1e6:.1f} Mb", fontsize=7, ha="right", va="top")
    ax.set_title(title, fontsize=9.5, fontweight="bold", pad=12)
cax = fig.add_axes([0.35, 0.05, 0.3, 0.018]); cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("mean genotype $r^2$ (binned); orange bar = block or inversion", fontsize=8)
fig.subplots_adjust(hspace=0.25, wspace=0.05, bottom=0.12)
fig.savefig(OUT, dpi=200, bbox_inches="tight"); print("wrote", OUT)
