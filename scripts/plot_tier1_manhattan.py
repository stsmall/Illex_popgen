"""Supp: genome-wide Manhattan of the diploSHIC-HMM (Tier-1) sweep signal. Per-100 kb
window empirical-outlier sweep score (S_pct) across all 45 autosomes + chr2/chr42; the
Tier-1 Markov calls are highlighted (hard vs soft), and the background-selection-robust
subset is marked. Shows where the diploSHIC-HMM calls fall genome-wide."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/scripts")
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/figS18_tier1_manhattan.png"
HARD, SOFT, ROBUST = "#D55E00", "#E69F00", "#000000"

win = pd.read_csv(f"{D}/outlier_scan_45_masked/windows.tsv", sep="\t")
win["chrom"] = win["chrom"].astype(str)
# append chrZ windows (male-based n=330 diploSHIC scan; diploSHIC-only, no RAiSD/SF2 on Z)
winZ = pd.read_csv(f"{D}/outlier_scan_Z/windows.tsv", sep="\t")
winZ["chrom"] = winZ["chrom"].astype(str)
win = pd.concat([win, winZ], ignore_index=True)
t1 = pd.read_csv(f"{D}/hmm_decode_masked/tier1_markov_calls.tsv", sep="\t")
t1["chrom"] = t1["chrom"].astype(str)
# append chrZ Tier-1 Markov calls (diploSHIC-only Tier-1; no BGS ring, no concordance)
t1Z = pd.read_csv(f"{D}/hmm_decode_masked/tier1_markov_calls_Z.tsv", sep="\t")
t1Z["chrom"] = t1Z["chrom"].astype(str)
t1 = pd.concat([t1, t1Z], ignore_index=True)
rob = pd.read_csv(f"{D}/hmm_decode_masked/tier1_BGS_robust_calls.tsv", sep="\t")
rob["chrom"] = rob["chrom"].astype(str)

# cumulative genome coordinate
chroms = sorted(win["chrom"].unique(), key=lambda c: int(c) if c.isdigit() else 99)
offset, cum, ticks, tlab = {}, 0, [], []
GAP = 5_000_000  # inter-chromosome gap so the short chrZ separates from chr45 and is not clipped
for c in chroms:
    offset[c] = cum
    L = win[win["chrom"] == c]["end"].max()
    ticks.append(cum + L / 2); tlab.append(c)
    cum += L + GAP
cum -= GAP
win["gpos"] = win.apply(lambda r: offset[r["chrom"]] + r["start"], axis=1)


def gp(df):
    return df.apply(lambda r: offset[r["chrom"]] + (r["start"] + r["end"]) / 2, axis=1)


fig, ax = plt.subplots(figsize=(14, 4.2))
# background windows, alternating grey shades by chromosome
for i, c in enumerate(chroms):
    w = win[win["chrom"] == c]
    ax.scatter(w["gpos"], w["S_pct"], s=3, color=("#C8C8C8" if i % 2 == 0 else "#9AA0A6"),
               linewidths=0, zorder=1, rasterized=True)
# Tier-1 calls (hard/soft), plotted at their peak S_pct
t1["gpos"] = gp(t1); t1["yval"] = t1["max_Spct"]
for ev, col, lab in [("hard", HARD, "Tier-1 hard"), ("soft", SOFT, "Tier-1 soft")]:
    s = t1[t1["evidence"] == ev]
    ax.scatter(s["gpos"], s["yval"], s=16, color=col, linewidths=0, zorder=3, label=f"{lab} ({len(s)})")
# BGS-robust subset: open ring at the call's own score (merge to get max_Spct)
robm = rob.merge(t1[["chrom", "start", "gpos", "max_Spct"]], on=["chrom", "start"], how="left")
robm = robm.dropna(subset=["gpos", "max_Spct"])
ax.scatter(robm["gpos"], robm["max_Spct"], s=46, facecolors="none", edgecolors=ROBUST,
           linewidths=0.9, zorder=4, label=f"BGS-robust ({len(robm)})")

ax.set_xticks(ticks); ax.set_xticklabels(tlab, fontsize=6.0)
ax.set_xlim(-cum * 0.005, cum * 1.005); ax.set_ylim(0.5, 1.005)
ax.set_xlabel("chromosome"); ax.set_ylabel("sweep score (empirical outlier percentile)")
ax.set_title("Genome-wide diploSHIC-HMM (Tier-1) sweep calls", loc="left", fontweight="bold", pad=8)
ax.legend(loc="upper center", ncol=3, fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, -0.13))
despine(ax)
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT, "| tier1:", len(t1), "robust:", len(rob))
