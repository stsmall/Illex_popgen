"""Separate GSEA-style dot plots for the diploSHIC-HMM (Tier-1) sweep-gene GO enrichment,
one per gene set: full Tier-1, background-selection-robust Tier-1, and hard-sweep Tier-1.
Each dot = an enriched biological-process term; x = fold enrichment, size = # genes,
colour = -log10(FDR). Top terms by FDR (the two themes surface naturally)."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs/hmm_decode"
SUP = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures"

SETS = [
    ("tier1_GO_enrichment.tsv",     "supp_gsea_tier1.png",     "Tier-1 (all diploSHIC-HMM calls)"),
    ("tier1_BGS_GO_enrichment.tsv", "supp_gsea_bgsrobust.png", "Tier-1, background-selection-robust"),
    ("tier1_GO_enrichment_hard.tsv","supp_gsea_hard.png",      "Tier-1 hard-sweep calls"),
]


def dotplot(infile, outfile, title, ntop=16):
    d = pd.read_csv(f"{D}/{infile}", sep="\t")
    ncol = "n_tier1" if "n_tier1" in d.columns else ("n_set" if "n_set" in d.columns else None)
    if ncol and ncol != "n_tier1":
        d = d.rename(columns={ncol: "n_tier1"})
    d = d[(d.FDR < 0.05) & (d.category == "biological_process")].copy()
    if len(d) == 0:
        print(f"skip {infile}: no significant terms"); return
    # de-duplicate near-identical parent terms by keeping the most significant, then top N by FDR
    d = d.sort_values("FDR").drop_duplicates("term").head(ntop).sort_values("fold")
    d["nlfdr"] = -np.log10(d.FDR)
    fig, ax = plt.subplots(figsize=(7.6, 0.34 * len(d) + 1.4))
    y = np.arange(len(d))
    smin, smax = d.n_tier1.min(), d.n_tier1.max()
    sizes = 30 + (d.n_tier1 - smin) / (smax - smin + 1e-9) * 240
    sc = ax.scatter(d.fold, y, s=sizes, c=d.nlfdr, cmap="viridis", edgecolors="white",
                    linewidths=0.6, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels([t[:36] for t in d.term], fontsize=8)
    ax.set_xlabel("fold enrichment (vs genome background)")
    ax.axvline(1.0, ls=":", lw=0.8, color=C["muted"], zorder=0)
    ax.set_xlim(1.0, d.fold.max() * 1.10)
    ax.set_title(f"GO enrichment — {title}", loc="left", fontweight="bold", fontsize=10.5)
    despine(ax)
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label(r"$-\log_{10}$ FDR", fontsize=8.5); cb.ax.tick_params(labelsize=7)
    for nn in [smin, int(np.median(d.n_tier1)), smax]:
        sz = 30 + (nn - smin) / (smax - smin + 1e-9) * 240
        ax.scatter([], [], s=sz, c="#888", edgecolors="white", linewidths=0.6, label=f"{int(nn)}")
    ax.legend(loc="lower right", fontsize=7.5, title="genes", title_fontsize=8,
              frameon=False, labelspacing=1.0, borderpad=0.6)
    fig.savefig(f"{SUP}/{outfile}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outfile}  ({len(d)} terms shown; {(pd.read_csv(f'{D}/{infile}', sep=chr(9)).FDR < 0.05).sum()} total FDR<0.05)")


for infile, outfile, title in SETS:
    dotplot(infile, outfile, title)
print("done")
