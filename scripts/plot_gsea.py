"""GO enrichment of the diploSHIC-HMM (Tier-1) sweep genes -- ONE consolidated figure.

Replaces the earlier set of near-identical GSEA dot plots (supp_gsea, supp_gsea_tier1,
supp_gsea_bgsrobust, supp_gsea_hard), which were confusing because "tier1" and "tier1"
were the same set and "hard"/"bgs" read like parallel categories when they are a
sub-classification and a screen.

This is a PAIRED dot plot on a single shared term axis:
  * blue circle    = all Tier-1 sweep genes (658 diploSHIC-HMM calls -> 476 GO-annotated)
  * orange diamond = the background-selection-robust subset (291 calls) -- the genes whose
                     local diversity reduction exceeds the background-selection expectation.
A grey arrow shows how each term's fold enrichment shifts between the two sets; the BGS
screen consistently SHARPENS the signal (rightward shift), and the light-response /
photoreception terms move the most. Marker area is proportional to the number of sweep
genes carrying the term. Only biological-process terms significant (FDR<0.05) in the full
Tier-1 set are shown; the hard-sweep subset is a sub-classification reported in the calls
table, not a separate enrichment.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs/hmm_decode"
OUT = os.path.join(_HERE, "supp_figures", "figS15_gsea.png")
COL_T1, COL_BG = C["accent"], C["warm"]   # blue = all Tier-1, orange = BGS-robust

# generic parent terms that carry no biology, plus deprecated GO ("obsolete ...") --
# drop so the two real themes (neural development; light response) read cleanly
DROP = {"regulation of cellular process", "cellular developmental process",
        "regulation of response to stimulus", "cellular anatomical entity morphogenesis",
        "regulation of developmental process", "regulation of biological process",
        "negative regulation of cellular process", "regulation of cellular process",
        "cellular process", "biological_process"}
NTOP = 15


def load(fn):
    d = pd.read_csv(f"{D}/{fn}", sep="\t")
    d = d[(d.category == "biological_process") & (d.FDR < 0.05)].copy()
    return d.set_index("term")


t1 = load("tier1_GO_enrichment.tsv")
bg = load("tier1_BGS_GO_enrichment.tsv")


def keep(t):
    return (t in t1.index and t in bg.index and t not in DROP
            and not t.startswith("obsolete "))


# rank terms significant in BOTH sets by their BEST FDR across the two sets, so the
# strongly BGS-sharpened light/photoreception terms surface alongside the neural ones
cand = [t for t in set(t1.index) & set(bg.index) if keep(t)]
best = {t: min(t1.loc[t, "FDR"], bg.loc[t, "FDR"]) for t in cand}
top = sorted(cand, key=lambda t: best[t])[:NTOP]
terms = sorted(top, key=lambda t: bg.loc[t, "fold"])   # ladder by BGS-robust fold

y = np.arange(len(terms))
f1 = t1.loc[terms, "fold"].to_numpy()
fb = bg.loc[terms, "fold"].to_numpy()
n1 = t1.loc[terms, "n_tier1"].to_numpy()
nb = bg.loc[terms, "n_tier1"].to_numpy()


def sz(n, nmin, nmax):
    return 40 + (n - nmin) / (nmax - nmin + 1e-9) * 300


nmn, nmx = min(n1.min(), nb.min()), max(n1.max(), nb.max())

fig, ax = plt.subplots(figsize=(8.4, 0.42 * len(terms) + 1.7))
# connector: Tier-1 -> BGS-robust (shows the sharpening)
for yi, a, b in zip(y, f1, fb):
    ax.annotate("", xy=(b, yi), xytext=(a, yi),
                arrowprops=dict(arrowstyle="-|>", color=C["muted"], lw=1.0,
                                shrinkA=6, shrinkB=6, alpha=0.8), zorder=1)
ax.scatter(f1, y, s=sz(n1, nmn, nmx), marker="o", facecolor=COL_T1, edgecolor="white",
           linewidths=0.7, zorder=3)
ax.scatter(fb, y, s=sz(nb, nmn, nmx), marker="D", facecolor=COL_BG, edgecolor="white",
           linewidths=0.7, zorder=4)

ax.axvline(1.0, ls=":", lw=0.9, color=C["muted"], zorder=0)
import textwrap
ax.set_yticks(y)
ax.set_yticklabels(["\n".join(textwrap.wrap(t, 30)) for t in terms], fontsize=8.5)
ax.set_ylim(-0.7, len(terms) - 0.3)
ax.set_xlim(0.9, max(f1.max(), fb.max()) * 1.12)
ax.set_xlabel("fold enrichment (vs. genome background)")
ax.set_title("GO enrichment of the diploSHIC-HMM sweep genes",
             loc="left", fontweight="bold", fontsize=11, pad=12)
despine(ax)

# --- legends: series identity + gene-count size ---
set_leg = [Line2D([], [], marker="o", ls="", mfc=COL_T1, mec="white", ms=9,
                  label="all Tier-1 (658 calls)"),
           Line2D([], [], marker="D", ls="", mfc=COL_BG, mec="white", ms=8,
                  label="BGS-robust subset (291)")]
l1 = ax.legend(handles=set_leg, loc="lower right", fontsize=8.5, frameon=False,
               title="sweep-gene set", title_fontsize=9, borderpad=0.6,
               bbox_to_anchor=(1.0, 0.02))
ax.add_artist(l1)
size_leg = [Line2D([], [], marker="o", ls="", mfc=C["muted"], mec="white",
                   ms=np.sqrt(sz(v, nmn, nmx)), label=f"{int(v)}")
            for v in (nmn, int(np.median([nmn, nmx])), nmx)]
ax.legend(handles=size_leg, loc="lower right", fontsize=8, frameon=False, title="# genes",
          title_fontsize=8.5, labelspacing=1.2, borderpad=0.7, bbox_to_anchor=(0.74, 0.02))

fig.savefig(OUT, dpi=220, bbox_inches="tight")
plt.close(fig)
print("wrote", OUT, f"| {len(terms)} paired terms | Tier-1 sig={len(t1)} BGS sig={len(bg)}")
