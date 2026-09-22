"""Supp Fig S23: the chr1 elevated-diversity block is shared by every sampling division and is not
differentiated among them. (a) pi per NAFO division (ANGSD, 50-kb windows) across chr1:19-36 Mb, lines
coloured by division latitude; (b) between-division Hudson FST (all 45 division pairs, ratio of averages)
in 100-kb windows over the same span; (c) block/flank pi ratio per division for all 22 blocks."""
import os, sys
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
_HERE = os.path.dirname(os.path.abspath(__file__)); _ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
B = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_perpop"
P = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/per_population/tables/perpop_windows_expanded.tsv"
M = "/sietch_colab/data_share/illex/popgen_data/seq_data/baker_2025/docs/Squid_Meta_Sept2024_Simple.txt"
OUT = os.path.join(_ROOT, "supp_figures", "figS23_pi_blocks_by_division.png")
S, E, LO, HI = 23.98, 31.28, 19.0, 36.3
lat = pd.read_csv(M, sep="\t").groupby("Zone").Lat.mean()
p = pd.read_csv(P, sep="\t"); p["chrom"] = p.chrom.astype(str)
p = p[(p.chrom == "1") & (p.WinCenter / 1e6 >= LO) & (p.WinCenter / 1e6 <= HI)]
f = pd.read_csv(f"{B}/block_fst_windows.tsv", sep="\t"); f = f[f.chrom.astype(str) == "1"]
r = pd.read_csv(f"{B}/block_pi_by_division.tsv", sep="\t")
divs = sorted(p.division.unique(), key=lambda d: lat[d])
norm = mpl.colors.Normalize(lat[divs].min(), lat[divs].max()); cm = mpl.cm.viridis

fig = plt.figure(figsize=(12, 8.2))
gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 0.8, 1.1], hspace=0.45)
ax = fig.add_subplot(gs[0])
ax.axvspan(S, E, color=C["shade"], lw=0, zorder=0)
for d in divs:
    q = p[p.division == d].sort_values("WinCenter")
    ax.plot(q.WinCenter / 1e6, q.pi.rolling(5, center=True, min_periods=1).mean(), lw=0.9, color=cm(norm(lat[d])), alpha=0.9)
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cm); cax = ax.inset_axes([0.73, 0.88, 0.22, 0.05])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.set_label("division mean latitude (°N)", fontsize=7); cb.ax.tick_params(labelsize=7)
ax.set_xlim(LO, HI); ax.set_ylabel(r"nucleotide diversity $\pi$")
ax.set_title("a   chr1 elevated-diversity block, each of ten sampling divisions separately", loc="left", fontweight="bold", fontsize=10)
despine(ax)
ax = fig.add_subplot(gs[1])
ax.axvspan(S, E, color=C["shade"], lw=0, zorder=0)
ax.axhline(0, color=C["muted"], lw=0.7, ls=":")
ax.plot((f.win + 5e4) / 1e6, f.fst, color=C["ink"], lw=0.8)
ins = f[f.region == "block"]; fl = f[f.region == "flank"]
ax.text(0.995, 0.92, f"between-division $F_{{ST}}$: block {ins.num.sum()/ins.den.sum():+.4f}, flanks {fl.num.sum()/fl.den.sum():+.4f}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8, color=C["ink"])
ax.set_xlim(LO, HI); ax.set_xlabel("chromosome 1 position (Mb)"); ax.set_ylabel(r"$F_{ST}$ (45 pairs)")
ax.set_title("b   no differentiation among divisions inside the block", loc="left", fontweight="bold", fontsize=10)
despine(ax)
ax = fig.add_subplot(gs[2])
order = r.groupby("block").ratio.median().sort_values().index.tolist()
for i, b in enumerate(order):
    q = r[r.block == b]
    ax.scatter(np.full(len(q), i) + np.linspace(-0.25, 0.25, len(q)), q.ratio, s=9,
               c=[cm(norm(lat[d])) for d in q.division], lw=0, zorder=3)
ax.axhline(1, color=C["muted"], lw=0.7, ls=":")
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=60, ha="right", fontsize=6.5)
ax.set_ylabel(r"$\pi$ block / $\pi$ flanks"); ax.set_ylim(0, None)
ax.set_title("c   all 22 blocks: the same excess in every division (points = divisions, coloured by latitude)", loc="left", fontweight="bold", fontsize=10)
despine(ax)
fig.savefig(OUT, dpi=200, bbox_inches="tight"); print("wrote", OUT)
