"""Fig: the recomb-aware Markov sweep decoder (design-informed HMM over the diploSHIC
per-window classification), + its two-tier agreement with the cross-method scan.

a. The model in action at a representative concordant sweep (chr6 ~12.7-13.1 Mb): the
   diploSHIC center posterior pC per 100 kb window, and the HMM Viterbi state track
   (neutral / linked / center) that localises the sweep centre inside the linked envelope.
b. Markov correction: raw diploSHIC hard/soft window calls -> spatially-coherent retained
   calls (isolated calls dropped), hard/soft split.
c. Two-tier agreement: fraction of windows decoded as sweep-centre at the Tier-2
   cross-method-concordant loci vs the genome background.
"""
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Patch

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
R = f"{D}/results/empirical_scan_fullsfs"
HMM = f"{R}/hmm_decode"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/fig_markov.png"

HARD = "#D55E00"   # vermillion
SOFT = "#E69F00"   # orange
LINK = "#BBBBBB"
NEUT = "#EDEDED"
STATE_COL = {"C": HARD, "LL": LINK, "LR": LINK, "N": "#8A8A8A", "GAP": "white"}

win = pd.read_csv(f"{R}/outlier_scan_45/windows.tsv", sep="\t")
win["chrom"] = win["chrom"].astype(str)
win["pC"] = win["p_hard"] + win["p_soft"]
sp = pd.read_csv(f"{HMM}/state_path.tsv.gz", sep="\t")
sp["chrom"] = sp["chrom"].astype(str)

# ---------- summary numbers (from tier1_markov.py) ----------
RAW_HARD, RAW_SOFT = 1094, 664           # raw diploSHIC argmax hard / soft windows
MK_HARD, MK_SOFT = 500, 217              # Markov-retained C windows, by peak pC
CONC_RETAIN, BG_RETAIN = 0.89, 0.029     # % windows Markov-C at concordant loci vs genome
OR_ENR, P_ENR = 282.0, 1.07e-47
N_CONC_REC, N_CONC = 30, 34

fig = plt.figure(figsize=(12.5, 4.2))
gs = gridspec.GridSpec(1, 3, width_ratios=[1.55, 1.0, 1.0], wspace=0.34,
                       left=0.06, right=0.985, top=0.88, bottom=0.14)

# ===== panel a: model in action, chr6 zoom (main bars + state strip) =====
gsa = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 0],
                                       height_ratios=[6, 0.7], hspace=0.08)
axa = fig.add_subplot(gsa[0])
axt = fig.add_subplot(gsa[1], sharex=axa)
c, lo, hi = "6", 11_000_000, 15_000_000
wc = win[(win.chrom == c) & (win.start >= lo) & (win.end <= hi)].sort_values("start")
sc = sp[(sp.chrom == c) & (sp.start >= lo) & (sp.start <= hi)].sort_values("start")
xm = (wc.start.to_numpy() + wc.end.to_numpy()) / 2 / 1e6
axa.bar(xm, wc.pC.to_numpy(), width=0.09, color=LINK, edgecolor="none", zorder=2)
stmap = dict(zip(sc.start, sc.state))
isC = wc.start.map(lambda s: stmap.get(s) == "C").to_numpy()
axa.bar(xm[isC], wc.pC.to_numpy()[isC], width=0.09,
        color=[HARD if h >= s else SOFT for h, s in zip(wc.p_hard[isC], wc.p_soft[isC])],
        edgecolor="none", zorder=3)
axa.set_ylim(0, 1.05); axa.set_ylabel(r"diploSHIC $p_{\rm sweep}$")
axa.set_title("a   Markov decoder at a concordant sweep", loc="left", fontweight="bold")
axa.tick_params(labelbottom=False)
despine(axa)
leg = [Patch(fc=HARD, label="centre: hard"), Patch(fc=SOFT, label="centre: soft"),
       Patch(fc=LINK, label="linked"), Patch(fc="#8A8A8A", label="neutral")]
axa.legend(handles=leg, fontsize=7.2, loc="upper right", ncol=2, handlelength=1.1,
           columnspacing=0.9, borderaxespad=0.3, frameon=False)
# state strip
for _, r in sc.iterrows():
    axt.add_patch(plt.Rectangle((r.start/1e6, 0), (r.end - r.start)/1e6, 1,
                                color=STATE_COL.get(r.state, "white"), lw=0))
axt.set_ylim(0, 1); axt.set_yticks([0.5]); axt.set_yticklabels(["HMM state"], fontsize=7.5)
axt.tick_params(axis="y", length=0)
axt.set_xlabel("chromosome 6 position (Mb)")
axt.set_xlim(lo/1e6, hi/1e6)
for s in axt.spines.values():
    s.set_visible(False)

# ===== panel b: Markov correction (stacked bars raw -> retained) =====
axb = fig.add_subplot(gs[0, 1])
x = [0, 1]
axb.bar(x, [RAW_HARD, MK_HARD], color=HARD, label="hard", zorder=2)
axb.bar(x, [RAW_SOFT, MK_SOFT], bottom=[RAW_HARD, MK_HARD], color=SOFT, label="soft", zorder=2)
for xi, tot in zip(x, [RAW_HARD+RAW_SOFT, MK_HARD+MK_SOFT]):
    axb.text(xi, tot + 30, f"{tot}", ha="center", va="bottom", fontsize=9, fontweight="bold")
axb.annotate("", xy=(0.78, MK_HARD+MK_SOFT+120), xytext=(0.22, RAW_HARD+RAW_SOFT+120),
             arrowprops=dict(arrowstyle="->", color=C["muted"], lw=1.3))
axb.text(0.5, RAW_HARD+RAW_SOFT+220, "−59% isolated", ha="center", fontsize=8, color=C["muted"])
axb.set_xticks(x); axb.set_xticklabels(["raw\ndiploSHIC", "Markov\nretained"], fontsize=8.5)
axb.set_ylabel("hard/soft window calls")
axb.set_ylim(0, 2200)
axb.set_title("b   Markov correction", loc="left", fontweight="bold")
axb.legend(fontsize=8, loc="upper right", frameon=False, handlelength=1.0)
despine(axb); axb.tick_params(axis="x", length=0)

# ===== panel c: concordance enrichment =====
axc = fig.add_subplot(gs[0, 2])
bars = axc.bar([0, 1], [CONC_RETAIN*100, BG_RETAIN*100], color=[HARD, LINK],
               width=0.62, zorder=2)
axc.set_xticks([0, 1])
axc.set_xticklabels(["concordant\nloci", "genome\nbackground"], fontsize=8.5)
axc.set_ylabel("windows decoded as sweep centre (%)")
axc.set_ylim(0, 100)
for b, v in zip(bars, [CONC_RETAIN*100, BG_RETAIN*100]):
    axc.text(b.get_x()+b.get_width()/2, v + 2, f"{v:.0f}%", ha="center", va="bottom",
             fontsize=9, fontweight="bold")
axc.set_title("c   Two-tier agreement", loc="left", fontweight="bold")
axc.text(0.5, 0.80, f"OR = {OR_ENR:.0f}\n$p<10^{{-46}}$\n{N_CONC_REC}/{N_CONC} concordant\nrecovered",
         transform=axc.transAxes, ha="center", va="center", fontsize=8, color=C["ink"],
         bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=C["muted"], lw=0.6))
despine(axc); axc.tick_params(axis="x", length=0)

fig.savefig(OUT, dpi=220, bbox_inches="tight")
print("wrote", OUT)
