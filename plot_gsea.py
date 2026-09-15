"""Supp: GSEA-style dot plot of the Tier-1 (diploSHIC-only, Markov-corrected) sweep-gene
GO enrichment. Each dot = one enriched biological-process term; x = fold enrichment,
dot size = number of Tier-1 genes, colour = -log10(FDR). Curated to the two dominant
themes (neural development; RNA-regulatory/transcriptional control)."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs/hmm_decode"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/supp_gsea.png"

d = pd.read_csv(f"{D}/tier1_GO_enrichment.tsv", sep="\t")
d = d[(d.FDR < 0.05) & (d.category == "biological_process")].copy()
# curated, ordered top -> shows both themes
want = ["generation of neurons", "neurogenesis", "neuron differentiation",
        "regulation of neurogenesis", "axonogenesis", "cell projection morphogenesis",
        "animal organ morphogenesis", "cell morphogenesis", "cell differentiation",
        "regulation of cell differentiation", "regulation of RNA metabolic process",
        "regulation of RNA biosynthetic process",
        "regulation of nucleic acid-templated transcription"]
rows = [d[d.term.str.lower() == t.lower()].iloc[0] for t in want if (d.term.str.lower() == t.lower()).any()]
s = pd.DataFrame(rows).drop_duplicates("go_id")
s = s.sort_values("fold")
s["nlfdr"] = -np.log10(s.FDR)

fig, ax = plt.subplots(figsize=(7.4, 4.6))
y = np.arange(len(s))
sizes = 30 + (s.n_tier1 - s.n_tier1.min()) / (s.n_tier1.max() - s.n_tier1.min() + 1e-9) * 240
sc = ax.scatter(s.fold, y, s=sizes, c=s.nlfdr, cmap="viridis", edgecolors="white",
                linewidths=0.6, zorder=3)
ax.set_yticks(y)
ax.set_yticklabels([t[:34] for t in s.term], fontsize=8)
ax.set_xlabel("fold enrichment (Tier-1 vs genome background)")
ax.axvline(1.0, ls=":", lw=0.8, color=C["muted"], zorder=0)
ax.set_xlim(1.0, s.fold.max() * 1.08)
ax.set_title("Tier-1 sweep genes: GO enrichment (biological process)", loc="left", fontweight="bold")
despine(ax)
cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
cb.set_label(r"$-\log_{10}$ FDR", fontsize=8.5); cb.ax.tick_params(labelsize=7)
# size legend
for nn, lab in [(s.n_tier1.min(), f"{int(s.n_tier1.min())}"), (int(s.n_tier1.median()), f"{int(s.n_tier1.median())}"), (s.n_tier1.max(), f"{int(s.n_tier1.max())} genes")]:
    sz = 30 + (nn - s.n_tier1.min()) / (s.n_tier1.max() - s.n_tier1.min() + 1e-9) * 240
    ax.scatter([], [], s=sz, c="#888", edgecolors="white", linewidths=0.6, label=lab)
ax.legend(loc="lower right", fontsize=7.5, title="Tier-1 genes", title_fontsize=8,
          frameon=False, labelspacing=1.0, borderpad=0.6)
fig.savefig(OUT, dpi=220, bbox_inches="tight")
print("wrote", OUT, "|", len(s), "terms")
