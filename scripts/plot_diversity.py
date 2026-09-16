"""Per-population diversity figure for the Illex illecebrosus manuscript.

Populations = NAFO sampling divisions. Everything is computed PER POPULATION
(not genome-wide-pooled).

MAIN (figures/fig6_diversity.png), three panels:
  a. Nucleotide diversity pi per population (box per division), ordered by latitude.
  b. Tajima's D per population (box per division), ordered by latitude.
  c. Pairwise Hudson FST among divisions -- clustermap (heatmap + hierarchical
     dendrogram). Near-zero everywhere -> no geographic structure (that is the result).

SUPP (supp_figures/figS9_diversity_boxplots.png):
  per-population, per-chromosome box plots of pi and Tajima's D.

Data:
  pi / Tajima's D : ROUTE A (ANGSD GL, folded SAF, accessible sites) -- correct for
                    this low-coverage unphased data (reproduces genome-wide pi~0.0093,
                    D~-2). steps/08_demography/per_population/tables/perpop_*.tsv
  FST matrix      : ROUTE B (pg_gpu Hudson FST on baker SNP callset), aggregated over
                    autosomes (chr2 inversion + chr42 sex + chrZ excluded). FST is a
                    variance ratio -> robust to the ascertainment bias. tables/fst_matrix.tsv
"""
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, SEQ_HUE, despine
apply()
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec, cm, colors
import matplotlib as mpl
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram

BASE = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/per_population/tables"
FIG = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/fig6_diversity.png"
SUP = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/figS9_diversity_boxplots.png"

# genome-wide ANGSD reference values (RESULTS_demography.md)
GW_PI, GW_D = 0.00930, -2.07

summ = pd.read_csv(f"{BASE}/perpop_diversity.tsv", sep="\t")
# prefer the expanded per-chromosome windows (more chromosomes) if present
import os as _os
_wf = f"{BASE}/perpop_windows_expanded.tsv" if _os.path.exists(f"{BASE}/perpop_windows_expanded.tsv") else f"{BASE}/perpop_windows.tsv"
win  = pd.read_csv(_wf, sep="\t")
print(f"windows source: {_os.path.basename(_wf)}")
fst  = pd.read_csv(f"{BASE}/fst_matrix.tsv", sep="\t", index_col=0)
win["chrom"] = win["chrom"].astype(str)

# latitude order (south -> north)
order = summ.sort_values("mean_lat")["division"].tolist()
lat = dict(zip(summ["division"], summ["mean_lat"]))

# latitude -> viridis colour (reinforces the ordering)
lmin, lmax = summ["mean_lat"].min(), summ["mean_lat"].max()
norm_lat = colors.Normalize(lmin, lmax)
cmap_lat = mpl.colormaps[SEQ_HUE]
lcol = {d: cmap_lat(norm_lat(lat[d])) for d in order}


BOX_FC = "#CBD5DC"   # single neutral fill (populations are labelled; no colour needed)
BOX_EC = "#5B6B75"


def boxpanel(ax, stat, ylabel, hline=None, ref=None, ref_label=None):
    data = [win[win["division"] == d][stat].to_numpy() for d in order]
    bp = ax.boxplot(data, positions=np.arange(len(order)), widths=0.66,
                    showfliers=False, patch_artist=True,
                    medianprops=dict(color=BOX_EC, lw=1.3))
    for patch in bp["boxes"]:
        patch.set(facecolor=BOX_FC, edgecolor=BOX_EC, lw=0.7)
    for w in ("whiskers", "caps"):
        for ln in bp[w]:
            ln.set(color=BOX_EC, lw=0.7)
    if hline is not None:
        ax.axhline(hline, color=C["muted"], ls="--", lw=0.8, zorder=0)
    if ref is not None:   # across-population median of the SAME windowed data (consistent basis)
        ax.axhline(ref, color=C["ink"], ls=":", lw=1.0, zorder=1)
        ax.text(len(order) - 0.45, ref, f"  {ref_label}", fontsize=7.5, color=C["ink"],
                va="center", ha="left")
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.6, len(order) - 0.4)
    despine(ax)
    ax.tick_params(axis="x", length=0)
    return bp


# ============ MAIN FIGURE ============
fig = plt.figure(figsize=(11.5, 5.0))
gs = gridspec.GridSpec(2, 2, width_ratios=[1.30, 1.0], height_ratios=[1, 1],
                       hspace=0.18, wspace=0.40,
                       left=0.085, right=0.955, top=0.9, bottom=0.11)

# across-population reference = median of ALL windowed values (same basis as the boxes,
# so it sits among them -- unlike the pooled genome-wide point estimate, which is deeper)
REF_PI = float(win["pi"].median())
REF_D = float(win["tajd"].median())

# --- panel a: pi ---
axa = fig.add_subplot(gs[0, 0])
boxpanel(axa, "pi", r"nucleotide diversity  $\pi$", ref=REF_PI, ref_label="median")
axa.set_xticks(np.arange(len(order)))
axa.tick_params(labelbottom=False)
axa.set_title("a", loc="left", fontweight="bold")

# --- panel b: Tajima's D ---
axb = fig.add_subplot(gs[1, 0], sharex=axa)
boxpanel(axb, "tajd", r"Tajima's  $D$", hline=0.0, ref=REF_D, ref_label="median")
axb.set_xticks(np.arange(len(order)))
axb.set_xticklabels([f"{d}\n{lat[d]:.0f}°N" for d in order], fontsize=8)
axb.set_title("b", loc="left", fontweight="bold")
axb.set_xlabel("NAFO division  (south → north)", labelpad=2)

# --- panel c: FST clustermap ---
pops = list(fst.index)
M = fst.loc[pops, pops].to_numpy().astype(float)
np.fill_diagonal(M, 0.0)
iu = np.triu_indices(len(pops), 1)
dist = M - M[iu].min()          # shift so smallest pairwise FST -> 0 (near-zero noise)
np.fill_diagonal(dist, 0.0)
Z = linkage(squareform(dist, checks=False), method="average")

gc = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[:, 1],
                                      height_ratios=[0.16, 1.0], width_ratios=[1.0, 0.045],
                                      hspace=0.03, wspace=0.05)
axd = fig.add_subplot(gc[0, 0])
dn = dendrogram(Z, labels=pops, ax=axd, color_threshold=0,
                above_threshold_color=C["muted"], no_labels=True)
leaves = dn["leaves"]
axd.set_axis_off()

leaf_pops = [pops[i] for i in leaves]
Mre = M[np.ix_(leaves, leaves)]
Mplot = Mre.copy(); np.fill_diagonal(Mplot, np.nan)   # blank diagonal

axh = fig.add_subplot(gc[1, 0])
im = axh.imshow(Mplot, cmap=SEQ_HUE, aspect="auto",
                vmin=np.nanmin(Mplot), vmax=np.nanmax(Mplot))
axh.set_xticks(np.arange(len(pops))); axh.set_yticks(np.arange(len(pops)))
axh.set_xticklabels(leaf_pops, fontsize=7.5, rotation=90)
axh.set_yticklabels(leaf_pops, fontsize=7.5)
axh.tick_params(length=0)
for s in axh.spines.values():
    s.set_visible(False)
axh.set_title("c", loc="left", fontweight="bold")

axcb = fig.add_subplot(gc[1, 1])
cbf = fig.colorbar(im, cax=axcb)
cbf.set_label(r"pairwise $F_{ST}$", fontsize=8.5)
cbf.ax.tick_params(labelsize=7)

fst_off = M[iu]
axh.text(0.5, -0.155,
         f"all $F_{{ST}}\\approx0$ ({fst_off.min():.3f} to {fst_off.max():.3f}); no geographic structure",
         transform=axh.transAxes, ha="center", va="top", fontsize=7.5, color=C["muted"])

fig.savefig(FIG, dpi=220, bbox_inches="tight")
print("wrote", FIG)
print("pi medians:", {d: round(float(win[win.division==d]['pi'].median()), 5) for d in order})
print("D  medians:", {d: round(float(win[win.division==d]['tajd'].median()), 3) for d in order})
print(f"FST off-diag: min={fst_off.min():.4f} max={fst_off.max():.4f} mean={fst_off.mean():.4f}")

# ============ SUPP FIGURE: per-population x per-chromosome ============
chroms = sorted(win["chrom"].unique(), key=lambda x: int(x) if str(x).isdigit() else 99)
nchrom = len(chroms)
chrom_cols = mpl.colormaps["cividis"](np.linspace(0.05, 0.9, nchrom))
figs, axs = plt.subplots(2, 1, figsize=(min(13, 3.0 + 1.05 * len(order)), 7.2), sharex=False)

for ax, stat, ylabel, hline, gw in ((axs[0], "pi", r"$\pi$", None, REF_PI),
                                    (axs[1], "tajd", "Tajima's $D$", 0.0, REF_D)):
    width = 0.82 / nchrom
    for ci, c in enumerate(chroms):
        data = [win[(win["division"] == d) & (win["chrom"] == c)][stat].to_numpy()
                for d in order]
        pos = np.arange(len(order)) + (ci - (nchrom - 1) / 2) * width
        bp = ax.boxplot(data, positions=pos, widths=width * 0.9, showfliers=False,
                        patch_artist=True, medianprops=dict(color="white", lw=0.6))
        for patch in bp["boxes"]:
            patch.set(facecolor=chrom_cols[ci], edgecolor="none", alpha=0.92)
        for w in ("whiskers", "caps"):
            for ln in bp[w]:
                ln.set(color=chrom_cols[ci], lw=0.5)
    if gw is not None:
        ax.axhline(gw, color=C["ink"], ls=":", lw=0.9, zorder=0)
    if hline is not None:
        ax.axhline(hline, color=C["muted"], ls="--", lw=0.8, zorder=0)
    ax.set_ylabel(ylabel); despine(ax); ax.tick_params(axis="x", length=0)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xlim(-0.6, len(order) - 0.4)

axs[0].set_xticklabels([])
axs[1].set_xticklabels([f"{d}\n{lat[d]:.0f}°N" for d in order], fontsize=8)
axs[1].set_xlabel("NAFO division  (south → north)")
handles = [plt.Rectangle((0, 0), 1, 1, color=chrom_cols[i]) for i in range(nchrom)]
axs[0].set_title("Per-population, per-chromosome diversity and Tajima's $D$", loc="left", pad=30)
axs[0].legend(handles, [f"chr{c}" for c in chroms], ncol=nchrom, fontsize=8,
              loc="lower center", bbox_to_anchor=(0.5, 1.02), frameon=False,
              handlelength=1.0, columnspacing=1.4)
figs.tight_layout()
figs.savefig(SUP, dpi=220, bbox_inches="tight")
print("wrote", SUP)
