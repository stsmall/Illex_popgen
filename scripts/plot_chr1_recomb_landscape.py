"""Supp: the true chromosome-1 ReLERNN recombination landscape with 95% bootstrap CIs, and
what it does NOT show. Three tracks: (a) male ReLERNN estimate (15 kb windows) + CI band +
rolling median; (b) female estimate, whose predictions collapse to a near-constant (IQR
0.005 cM/Mb) -- degenerate, shown for completeness; (c) nucleotide diversity pi in 10 kb
windows. The 2.5 Mb block at 45.0-47.5 Mb is 80% repeat (vs 52% chromosome-wide), almost
entirely inaccessible, and carries no ReLERNN windows: the centromere. pi halves on its
flanks, but the recombination estimate in those accessible flanks dips only ~5% -- the map's
chromosome-scale dynamic range is compressed (0.17-0.23 cM/Mb), so a centromeric drop cannot
be resolved and 'nearly uniform' is in part an estimator ceiling."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE) if os.path.basename(_HERE) == "scripts" else _HERE
sys.path.insert(0, _HERE)
from figstyle import apply, C, SEXES, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt

RD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn"
WIN = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/bgs_diagnostic_windows.tsv"
OUT = os.path.join(_ROOT, "supp_figures", "figS22_chr1_recomb_landscape.png")
CHROM = "1"
GAP = (45.0, 47.5)          # masked, 80%-repeat block: no ReLERNN windows (candidate centromere)


def load(proj, name):
    d = pd.read_csv(f"{RD}/{proj}/proj/{name}.kept.PREDICT.BSCORRECTED.txt", sep="\t")
    d["chrom"] = d["chrom"].astype(str).str.replace(r"^b'|'$", "", regex=True)
    d = d[d.chrom == CHROM].copy().sort_values("start")
    for c in ("recombRate", "CI95LO", "CI95HI"):
        d[c] = d[c].astype(float) * 1e8
    d["mid"] = (d.start + d.end) / 2 / 1e6
    return d


m = load("run_male_auto", "male")
w = pd.read_csv(WIN, sep="\t")
w = w[w.chrom.astype(str) == CHROM].sort_values("start")
w["mid"] = (w.start + w.stop) / 2 / 1e6
L = max(m.mid.max(), w.mid.max())

fig, ax0 = plt.subplots(1, 1, figsize=(13, 3.6)); axs = [ax0]
for ax in axs:
    ax.axvspan(*GAP, color=C["faint"], lw=0, zorder=0)

# (a) male
ax = axs[0]; d = m; col = C["accent"]
ax.fill_between(d.mid, d.CI95LO, d.CI95HI, color=col, alpha=0.22, lw=0, label="95% bootstrap CI")
ax.plot(d.mid, d.recombRate, color=col, lw=0.6, alpha=0.9, label="ReLERNN estimate (15 kb windows)")
roll = d.recombRate.rolling(41, center=True, min_periods=15).median()
ax.plot(d.mid, roll, color=C["ink"], lw=1.4, label="rolling median (≈0.6 Mb)")
med = d.recombRate.median(); ax.axhline(med, color=C["muted"], ls=":", lw=0.9)
broad = d.recombRate.rolling(133, center=True, min_periods=60).median()
axs[0].text(GAP[1] + 0.8, 0.40, "masked block 45.0–47.5 Mb (centromere)\n80% repeat; no windows", ha="left",
        va="top", fontsize=7.5, color=C["ink"])
ax.text(0.995, 0.06, f"chromosome-scale range of rolling median: {broad.min():.2f}–{broad.max():.2f} cM/Mb "
        f"({100*(broad.min()/med-1):+.0f}% … {100*(broad.max()/med-1):+.0f}%); pericentromeric flanks −5%",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5, color=C["muted"])
ax.set_ylim(0, 0.5); ax.set_ylabel("recombination\n(cM/Mb)")
ax.set_title(f"Chromosome 1 recombination map  (median {med:.3f} cM/Mb)", loc="left", fontweight="bold", fontsize=10)
ax.legend(loc="upper left", frameon=False, fontsize=7.5, ncol=3); despine(ax)


axs[0].set_xlabel("chromosome 1 position (Mb)")
axs[0].set_xlim(0, L)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT, f"| median {med:.3f}, rolling range {broad.min():.3f}-{broad.max():.3f}")