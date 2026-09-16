"""Fig 5A — selection-scan pipeline schematic (the 'selection Markov model' as a
workflow): simulate -> classify -> HMM-decode -> empirical-outlier + concordance
-> background-selection filter -> candidates."""
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C
apply()
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/fig_pipeline.png"

fig, ax = plt.subplots(figsize=(8.4, 4.6))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, title, sub, fc="#eef3f7", ec=C["accent"]):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12",
                 fc=fc, ec=ec, lw=1.3, zorder=2))
    ax.text(x + w / 2, y + h * 0.63, title, ha="center", va="center", fontsize=9.5,
            fontweight="bold", color=C["ink"], zorder=3)
    if sub:
        ax.text(x + w / 2, y + h * 0.26, sub, ha="center", va="center", fontsize=7.4,
                color="#555", zorder=3)

def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=13,
                 lw=1.4, color=C["muted"], zorder=1, shrinkA=2, shrinkB=2))

# --- top row: build the classifier ---
box(0.2, 7.4, 3.0, 2.1, "Simulate", "SLiM + msprime under\nfitted growth demography;\nfull SFS, real missingness",
    fc="#f3eef7", ec="#7d5ba6")
box(3.7, 7.4, 3.0, 2.1, "Train classifier", "diploSHIC 5-class CNN\n(neutral / hard / soft /\nlinked-hard / linked-soft)",
    fc="#f3eef7", ec="#7d5ba6")
box(7.2, 7.4, 2.6, 2.1, "Scan genome", "per-window class probs;\nsweep score\nS = p(hard)+p(soft)",
    fc="#f3eef7", ec="#7d5ba6")
arrow(3.2, 8.45, 3.7, 8.45); arrow(6.7, 8.45, 7.2, 8.45)

# down to filtering
arrow(8.5, 7.4, 8.5, 6.3)

# --- middle: the problem + Markov decode ---
box(5.6, 4.2, 4.2, 2.0, "Markov / HMM decode",
    "S over-calls (sim≠real spatial\nheterogeneity); decode the window\nchain, don't threshold raw S",
    fc="#fff4e8", ec=C["warm"])
arrow(7.7, 6.3, 7.7, 6.2)

# --- filters (left) ---
box(0.2, 4.2, 4.8, 2.0, "Empirical-outlier ranking",
    "top-percentile WITHIN\ncallability × gene-density strata\n(controls accessibility + BGS)",
    fc="#eef7f1", ec="#009E73")
arrow(5.6, 5.2, 5.0, 5.2)

box(0.2, 1.6, 4.8, 2.0, "Cross-method concordance",
    "require RAiSD or SweepFinder2\nto agree  →  discounts the\nclassifier over-call",
    fc="#eef7f1", ec="#009E73")
arrow(2.6, 4.2, 2.6, 3.6)

box(5.6, 1.6, 4.2, 2.0, "Background-selection test",
    "keep only π reduced beyond the\nlocal coding-density expectation",
    fc="#eef7f1", ec="#009E73")
arrow(5.0, 2.6, 5.6, 2.6)

# --- result ---
box(3.0, -0.3, 4.0, 1.2, "Candidate sweeps",
    "34 concordant  →  11 BGS-robust", fc="#fdecea", ec=C["warm"])
arrow(6.6, 1.6, 6.0, 0.9); arrow(3.4, 1.6, 4.0, 0.9)

fig.savefig(OUT, dpi=220, bbox_inches="tight")
print("wrote", OUT)
