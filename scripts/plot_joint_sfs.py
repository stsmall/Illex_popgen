"""Supp Fig S5: joint 2D SFS between the chr2 arrangements, inversion vs collinear control.

Input: unfolded realSFS 2D spectra (81 x 81, 40 2a/2a and 40 2+/2+ individuals),
m[i, j] = sites with i non-reference alleles in 2a and j in 2+, plus the summary
in sfs2d_karyotype.json (msinv repo, illex/scripts/sfs2d_karyotype.py).
Frequencies are of the NON-REFERENCE allele, not the derived allele: the
reference assembly carries 2+ in the inversion.
"""
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R = Path("/sietch_colab/ssmall/projects/msinv_dir/msinv_pboot_wt/results/illex")
OUT = Path(__file__).resolve().parent.parent / "supp_figures" / "figS5_joint_sfs.png"

stats = json.load(open(R / "sfs2d_karyotype.json"))
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
fig.subplots_adjust(wspace=0.32)
cmap = mpl.colormaps["magma"].copy()
cmap.set_bad("#d9d9d9")

for ax, (key, tag, label) in zip(axes, (("body", "a", "inversion  (60.54–79.50 Mb)"),
                                        ("control", "b", "collinear control  (10–30 Mb)"))):
    m = np.loadtxt(R / "sfs2d" / f"AA_BB.{key}.2dSFS").reshape(81, 81)
    v = m.copy()
    v[0, 0] = v[-1, -1] = np.nan                       # monomorphic corners
    lv = np.where(v > 0, np.log10(np.maximum(v, 1.0)), np.nan)
    im = ax.imshow(lv, origin="lower", cmap=cmap, aspect="equal", vmin=0, vmax=4.2,
                   extent=(0, 1, 0, 1), interpolation="bilinear")
    ax.plot([0, 1], [0, 1], lw=0.6, ls=(0, (2, 2)), color="#bdbdbd", alpha=0.7)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("non-reference allele frequency in 2+  (ancestral)")
    if key == "body":
        ax.set_ylabel("non-reference allele frequency in 2a  (derived)")
    s = stats[key]
    nf = s["nearfix_AAalt"] + s["nearfix_BBalt"]
    ax.set_title(f"{tag}   {label}", loc="left", fontweight="bold")
    ax.text(0.035, 0.965, f"shared polymorphism   {100 * s['shared'] / s['variable']:.1f}%\n"
            f"near-fixed differences   {nf:,.0f}", transform=ax.transAxes, ha="left", va="top",
            fontsize=8, color="white", linespacing=1.45,
            bbox=dict(fc="#18181b", ec="none", alpha=0.86, pad=3.2))

cb = fig.colorbar(im, ax=axes, fraction=0.030, pad=0.020)
cb.set_label(r"sites per cell  ($\log_{10}$, floored at 1)", fontsize=8)
cb.outline.set_visible(False)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
