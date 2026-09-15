"""Shared matplotlib style for all Illex manuscript figures.

Import at the top of every plotting script:
    import sys; sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
    from figstyle import apply, C, KARYO, SEXES, despine, save

One consistent look: Okabe-Ito colourblind-safe categorical palette (fixed
order, never cycled), thin recessive axes, minimal chartjunk, publication dpi.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito (CVD-safe, the genomics standard). Fixed order; do not cycle.
OKABE = ["#000000", "#E69F00", "#56B4E9", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

# semantic assignments (entity -> colour, never by rank)
C = {
    "ink":    "#222222",   # primary marks / text
    "muted":  "#8a8a8a",   # secondary / grid
    "faint":  "#d9d9d9",   # backgrounds / bands
    "accent": "#0072B2",   # single-series highlight (blue)
    "warm":   "#D55E00",   # emphasis / selection (vermillion)
    "shade":  "#f4a58244", # inversion region shading (translucent)
}
KARYO = {"AA": "#0072B2", "AB": "#E69F00", "BB": "#009E73"}   # ancestral/het/derived (blue/orange/green, CVD-safe)
SEXES = {"male": "#0072B2", "female": "#D55E00"}
SEQ_HUE = "viridis"        # sequential magnitude (perceptually uniform)


def apply():
    plt.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 220,
        "savefig.bbox": "tight", "savefig.facecolor": "white",
        "figure.facecolor": "white", "axes.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "legend.fontsize": 9, "legend.frameon": False,
        "axes.linewidth": 0.8, "axes.edgecolor": "#444444",
        "axes.grid": False, "grid.color": "#e6e6e6", "grid.linewidth": 0.6,
        "lines.linewidth": 1.4, "lines.solid_capstyle": "round",
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "axes.titlelocation": "left", "axes.titlepad": 8,
    })


def despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)


def save(fig, path):
    fig.savefig(path)
    print("wrote", path)
