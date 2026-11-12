#!/usr/bin/env python
"""Fig 2 -- the chr2 inversion as a genomic object (3-panel composite).

    /home/ssmall/miniforge3/envs/bioinfo-buddy/bin/python fig2_inversion.py

Panels
  a  chr2 karyotype PCA ("three stripes") -- PC1 vs PC2 coloured by karyotype,
     recovered as the true AA/AB/BB assignment via a 3-component Gaussian
     mixture on PC1 (ordered by component mean: low->AA, mid->AB, high->BB),
     fitted on all 653 PCA rows and then restricted to the 621 analysed individuals.
     Data: steps/03_karyotype/chr2_karyo_coords.tsv (+ sample IDs from karyotypes.tsv)
  b  F_ST across the whole of chr2 (AA vs BB), inversion shaded -- shows the
     inversion (60.54-79.50 Mb) is the only differentiated region.
     Data: steps/03_karyotype/inversion_content/chr2_karyo_fst_windows.tsv, the same 100-kb
     Hudson ratio-of-averages windows (all arrangement homozygotes) and the same filter
     (>=10 SNPs, denominator > 0, F_ST clipped at 0) as Supplementary Fig. S20
     (plot_chr2_karyotype_fst.py), so the two figures show identical data.
  c  per-arrangement diversity (d_xy, pi_BB ancestral, pi_AA derived) across
     the inversion plus 20 Mb of collinear sequence on each side (40-100 Mb) --
     shows the derived AA arrangement is the LESS diverse one and that the three
     curves collapse together outside the breakpoints.
     Data: steps/03_karyotype/inversion_content/flank_div/chr2_40_100_karyo_windows_snps.csv
     (compute_flank_div.py: pg_gpu windowed pi/dxy, 100 kb, AA n=254 / BB n=95)

Panel c reproduces the data + logic of ``_fst_chr2``/``_diversity`` and the
manu_fig_inv_fst / manu_fig_inv_diversity functions in
``/sietch_colab/ssmall/projects/msinv_dir/inversion_sims/files/illex/scripts/figures_inversion.py``,
re-plotted natively in the shared manuscript style rather than reusing the PNGs.
"""
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/scripts")
from figstyle import apply, C, KARYO, despine
apply()
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import os
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

KARYO_COORDS = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/chr2_karyo_coords.tsv"
# chr2_karyo_coords.tsv carries no sample IDs. Its rows are in the same order as karyotypes.tsv
# (identical PC1/PC2), which supplies them. Panel a shows the 621 analysed Baker et al. 2025
# individuals: the 20 non-Baker rows (population A/B/C) and the 12 BAMs with essentially no data
# (data/qc/empty_bams_12.tsv, which sit at the PCA origin) are dropped.
KARYO_IDS = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/karyotypes.tsv"
BAKER_633 = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/baker_633.txt"
EMPTY_BAMS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "qc", "excluded_14.tsv")
DISP = {"AA": "2a/2a", "AB": "2a/2+", "BB": "2+/2+"}   # display names (2a derived, 2+ ancestral/reference)
FST_WINDOWS = ("/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/"
               "chr2_karyo_fst_windows.tsv")   # same table as Supplementary Fig. S20
DIV_WINDOWS = ("/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/"
               "flank_div/chr2_40_100_karyo_windows_snps.csv")   # chr2:40-100 Mb, biallelic SNPs, same pg_gpu recipe as msinv
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures", "fig3_inversion.png")

BP_L, BP_R = 60_540_000, 79_500_000     # pinned inversion breakpoints (Mb)
CHR2_LEN = 119_466_599


# --------------------------------------------------------------------------
# panel a -- karyotype PCA (GMM on PC1 recovers the true AA/AB/BB call)
# --------------------------------------------------------------------------
def _karyotype_pca():
    d = pd.read_csv(KARYO_COORDS, sep="\t")
    ids = pd.read_csv(KARYO_IDS, sep="\t")
    assert len(ids) == len(d) and np.allclose(ids["PC1"], d["PC1"]) and np.allclose(ids["PC2"], d["PC2"])
    d["sample"] = ids["sample"].values
    # GMM fitted on all 653 rows of the PCA, as the karyotype calls were made (on this set it
    # reproduces karyotypes.tsv exactly). Refitting on the 621 alone does not, so the subset
    # is taken after labelling.
    x = d["PC1"].values.reshape(-1, 1)
    gmm = GaussianMixture(n_components=3, random_state=0, n_init=5).fit(x)
    labels = gmm.predict(x)
    order = np.argsort(gmm.means_.ravel())          # low -> high PC1 mean
    name_of = {order[0]: "AA", order[1]: "AB", order[2]: "BB"}
    d["karyotype"] = [name_of[l] for l in labels]
    assert (d["karyotype"].values == ids["karyotype"].values).all()
    baker = set(open(BAKER_633).read().split())
    empty = set(pd.read_csv(EMPTY_BAMS, sep="\t", comment="#")["sample"])
    d = d[d["sample"].isin(baker) & ~d["sample"].isin(empty)].reset_index(drop=True)
    print("panel a: n =", len(d), d["karyotype"].value_counts().reindex(["AA", "AB", "BB"]).to_dict())
    return d


# --------------------------------------------------------------------------
# panel b -- FST across chr2 (AA vs BB), 100 kb windows
# --------------------------------------------------------------------------
def _fst_chr2():
    # identical recipe to plot_chr2_karyotype_fst.py (Supplementary Fig. S20)
    w = pd.read_csv(FST_WINDOWS, sep="\t", header=None, names=["win", "num", "den", "n"])
    w = w[(w.n >= 10) & (w.den > 0)].sort_values("win")
    fst = (w.num / w.den).clip(lower=0).to_numpy()
    mid = (w.win.to_numpy() + 50_000)
    inv = (mid >= BP_L) & (mid <= BP_R)
    print(f"panel b: {len(w)} windows | median FST inversion={np.median(fst[inv]):.4f} "
          f"collinear={np.median(fst[~inv]):.4f} | inversion max={fst[inv].max():.4f}")
    return mid / 1e6, fst                              # window midpoint, Mb


# --------------------------------------------------------------------------
# panel c -- per-arrangement diversity across the inversion, 100 kb windows
# --------------------------------------------------------------------------
def _diversity():
    d = pd.read_csv(DIV_WINDOWS)
    d = d.sort_values("window_start")
    mid = d["window_start"].values / 1e6 + 0.05
    return mid, d["pi_AA"].values, d["pi_BB"].values, d["dxy"].values


def main():
    karyo = _karyotype_pca()
    fw, ff = _fst_chr2()
    mid, pAA, pBB, dxy = _diversity()

    fig = plt.figure(figsize=(9.0, 6.8))
    gs = GridSpec(2, 2, width_ratios=[1.0, 1.55], height_ratios=[1.0, 0.85],
                  hspace=0.42, wspace=0.32, figure=fig)

    # ---- a. karyotype PCA -------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    for k in ("AA", "AB", "BB"):
        d = karyo[karyo["karyotype"] == k]
        ax.scatter(d["PC1"], d["PC2"], s=10, color=KARYO[k], label=DISP[k],
                   alpha=0.85, linewidths=0)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_box_aspect(1)
    ax.set_title("a", loc="left", fontweight="bold")
    ax.legend(loc="upper center", ncol=3, handletextpad=0.3, borderaxespad=0.2,
              columnspacing=1.0, markerscale=1.5, frameon=True, facecolor="white",
              edgecolor="none", framealpha=0.85)
    despine(ax)

    # ---- b. FST across chr2 -----------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.axvspan(BP_L / 1e6, BP_R / 1e6, color=C["shade"], lw=0)
    ax.plot(fw, ff, lw=0.5, color=C["ink"])
    ax.set_xlim(0, CHR2_LEN / 1e6)
    ax.set_ylim(0, 0.65)
    ax.set_yticks([0, 0.2, 0.4, 0.6])
    ax.set_xlabel("chr2 position (Mb)")
    ax.set_ylabel(r"$F_{ST}$  (2a vs 2+)")
    ax.set_title("b", loc="left", fontweight="bold")
    despine(ax)

    # ---- c. per-arrangement diversity across the inversion ----------------
    ax = fig.add_subplot(gs[1, :])
    ax.axvspan(BP_L / 1e6, BP_R / 1e6, color=C["shade"], lw=0, zorder=0)
    COLINEAR_DXY = 1.270   # median d_xy(AA,BB) in matched colinear control windows (x1e-3)
    ax.axhline(COLINEAR_DXY, ls="--", lw=1.0, color=C["warm"], alpha=0.65, zorder=1,
               label=r"collinear $d_{xy}$ (matched controls, 10–30 Mb)")
    ax.plot(mid, dxy * 1e3, lw=1.1, color=C["warm"], label=r"$d_{xy}$ (2a, 2+)")
    ax.plot(mid, pBB * 1e3, lw=1.1, color=KARYO["BB"], label=r"$\pi_{2+}$ (ancestral)")
    ax.plot(mid, pAA * 1e3, lw=1.1, color=KARYO["AA"], label=r"$\pi_{2a}$ (derived)")
    ax.set_xlim(40, 100)
    ax.set_xticks([40, 50, 60, 70, 80, 90, 100])
    ymax = float(np.nanpercentile(dxy * 1e3, 99.8)) * 1.05
    ax.set_ylim(0, ymax)
    ax.set_xlabel("chr2 position (Mb)")
    ax.set_ylabel(r"diversity  ($\times 10^{-3}$)")
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=4, borderaxespad=0.0,
              columnspacing=1.4, handlelength=1.6, handletextpad=0.4, frameon=False)
    ax.set_title("c", loc="left", fontweight="bold", pad=22)
    despine(ax)

    fig.savefig(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
