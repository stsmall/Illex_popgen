import sys; sys.path.insert(0, '/sietch_colab/data_share/illex/popgen_data/analysis/manuscript')
from figstyle import apply, C, despine, save
apply()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import csv

# 45-chromosome outlier scan (includes chr2 and chr42; the older
# outlier_scan/windows.tsv is a stale 43-chromosome run missing both)
WINDOWS = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs/outlier_scan_45/windows.tsv"
# chrZ (male-based n=330 diploSHIC scan). RAiSD/SweepFinder2 do NOT cover chrZ, so
# Z windows are diploSHIC-only (n_methods<=1, never a >=2-method "confirmed candidate");
# the Z Tier-1 hard/soft calls are overlaid separately below.
WINDOWS_Z = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs/outlier_scan_Z/windows.tsv"
TIER1_Z = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/empirical_scan_fullsfs/hmm_decode/tier1_markov_calls_Z.tsv"
# chrZ Tier-1 calls are one category (diploSHIC-only, male-based) -> a single distinct
# colour (blue, the male convention used elsewhere), not the orange sweep-candidate family
HARD_Z = SOFT_Z = "#0072B2"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/fig7_selection_manhattan.png"

# ---- why the y-axis is NOT raw diploSHIC S ----
# Raw S = P(hard)+P(soft) saturates at ~1 across a huge, spatially-correlated
# fraction of the genome (a documented sim-vs-real mismatch: the CNN
# over-calls on real heterogeneous data). Plotting raw S makes thousands of
# non-candidate windows look identical to the true candidates, which is
# exactly the "why did so few pass" confusion this figure must avoid.
#
# Instead we plot each window's STRATIFIED rank: S_pct is diploSHIC S's
# percentile rank *within its own callability x gene-density stratum*
# (windows.tsv, columns S_pct/S_out) -- this is the value actually used to
# call outliers, and it discounts strata where S=1 is the common background
# rather than a real signal. We display it as -log10(1 - S_pct) purely so
# the top of the rank distribution (where the real outliers live) is spread
# out and reads as a Manhattan-style peak instead of a saturated ceiling.
S_OUT_PCT = 0.99  # windows.tsv's own stratified-outlier cut (top 1%/stratum)
EPS = 1e-6

def rank_y(pct):
    return -np.log10(max(1.0 - pct, EPS))

Y_THRESH = rank_y(S_OUT_PCT)

# ---- load per-window scan (autosomes + chrZ) ----
rows = []
for wfile in (WINDOWS, WINDOWS_Z):
    with open(wfile) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            chrom = row["chrom"]
            start = int(row["start"])
            end = int(row["end"])
            S_pct = float(row["S_pct"])
            n_methods = int(row["n_methods"]) if row["n_methods"] not in ("", "NA") else 0
            S_out = row["S_out"] == "1"
            candidate = S_out and n_methods >= 2   # stratified outlier AND
                                                    # corroborated by >=1 independent
                                                    # method (SF2/RAiSD) AND BGS-robust
                                                    # (candidate set already filtered
                                                    # upstream to n_methods>=2)
            rows.append(dict(chrom=chrom, start=start, end=end, S_pct=S_pct,
                              y=rank_y(S_pct), S_out=S_out, candidate=candidate))

def chrom_sort_key(c):
    return (0, int(c)) if c.isdigit() else (1, 99)

chroms = sorted({r["chrom"] for r in rows}, key=chrom_sort_key)
chrom_len = {c: max(r["end"] for r in rows if r["chrom"] == c) for c in chroms}

offset = {}
cum = 0
GAP = 5_000_000  # inter-chromosome gap so the short chrZ separates from chr45 and is not clipped
for c in chroms:
    offset[c] = cum
    cum += chrom_len[c] + GAP
cum -= GAP  # no trailing gap past the last chromosome

def gpos(chrom, pos):
    return offset[chrom] + pos

background = [r for r in rows if not r["S_out"]]
outlier_only = [r for r in rows if r["S_out"] and not r["candidate"]]
candidates = [r for r in rows if r["candidate"]]

print(f"genome-wide windows: {len(rows)} across {len(chroms)} chroms")
print(f"stratified outliers (top {100*(1-S_OUT_PCT):.0f}% per stratum): "
      f"{len(outlier_only) + len(candidates)}")
print(f"  -> corroborated (>=2 methods) + BGS-robust candidates: {len(candidates)}")

# ---- named genes to label ----
genes = [
    ("1", 20_300_000, "Aplnr"),
    ("35", 39_000_000, "ZEB2"),
    ("32", 9_300_000, "MACROD2"),
    ("5", 38_000_000, "cept1"),
    ("16", 35_300_000, "novel"),
]

fig, ax = plt.subplots(figsize=(9.5, 3.4))

# alternating chromosome bands
for i, c in enumerate(chroms):
    if i % 2 == 0:
        continue
    ax.axvspan(offset[c], offset[c] + chrom_len[c], color=C["faint"], lw=0, zorder=0)

# stratified-outlier reference line (drawn under the points); labelled by the y-axis
# 0.99 tick and the cascade note above the plot, so no floating label needed here
ax.axhline(Y_THRESH, color=C["muted"], lw=0.8, ls="--", zorder=1)

def xy(recs):
    x = [gpos(r["chrom"], (r["start"] + r["end"]) / 2) for r in recs]
    y = [r["y"] for r in recs]
    return x, y

# 1) discounted background -- the raw over-call lives in here, but ranked
#    within its own stratum it sits low, near everything else
bx, by = xy(background)
ax.scatter(bx, by, s=1.6, color=C["ink"], alpha=0.22, linewidths=0,
           zorder=2, rasterized=True)

# 2) stratified outliers that failed corroboration/BGS -- real rank outliers,
#    but not confirmed candidates
ox, oy = xy(outlier_only)
ax.scatter(ox, oy, s=7, color=C["muted"], alpha=0.75, linewidths=0,
           zorder=3, rasterized=True)

# 3) the corroborated, BGS-robust candidates
cx, cy = xy(candidates)
ax.scatter(cx, cy, s=24, color=C["warm"], edgecolors="white", linewidths=0.4,
           zorder=4)

# 3b) chrZ Tier-1 (diploSHIC-only) calls -- RAiSD/SF2 do not cover chrZ, so these
#     cannot enter the >=2-method concordant set; shown in Tier-1 hard/soft styling
#     (diamond markers) with no BGS ring.
zt1 = []
with open(TIER1_Z) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        zt1.append((gpos("Z", (int(row["start"]) + int(row["end"])) / 2),
                    rank_y(float(row["max_Spct"])), row["evidence"]))
if zt1:
    zx_h = [x for x, y, e in zt1 if e == "hard"]; zy_h = [y for x, y, e in zt1 if e == "hard"]
    zx_s = [x for x, y, e in zt1 if e == "soft"]; zy_s = [y for x, y, e in zt1 if e == "soft"]
    if zx_h:
        ax.scatter(zx_h, zy_h, s=26, marker="D", color=HARD_Z, edgecolors="white",
                   linewidths=0.4, zorder=5)
    if zx_s:
        ax.scatter(zx_s, zy_s, s=26, marker="D", color=SOFT_Z, edgecolors="white",
                   linewidths=0.4, zorder=5)

# label a handful of the strongest named candidate genes
label_offsets = [(0, 6), (0, 6), (-16, 6), (0, 6), (12, 6)]
for (chrom, pos, name), (dx, dy) in zip(genes, label_offsets):
    win_start = pos + 1
    hit = [r for r in rows if r["chrom"] == chrom and r["start"] == win_start]
    x = gpos(chrom, pos + 50_000)
    y = hit[0]["y"] if hit else Y_THRESH
    ax.annotate(name, xy=(x, y), xytext=(dx, dy), textcoords="offset points",
                ha="center", va="bottom", fontsize=7.5, style="italic",
                color=C["ink"], zorder=5)

# note on the filtering cascade -- placed ABOVE the plot (top margin) so it never covers points
n_out = len(outlier_only) + len(candidates)
ax.text(0.0, 1.20,
         f"{n_out} stratified outliers  →  {len(candidates)} corroborated windows in 34 regions "
         "(diploSHIC-HMM outlier + a footprint method; Tier-2)",
         transform=ax.transAxes, ha="left", va="bottom", fontsize=6.8,
         color=C["ink"])

# compact legend
legend_handles = [
    Line2D([0], [0], marker="o", linestyle="", markersize=3.5,
           markerfacecolor=C["ink"], markeredgewidth=0, alpha=0.5,
           label="discounted (in-stratum rank low)"),
    Line2D([0], [0], marker="o", linestyle="", markersize=4.5,
           markerfacecolor=C["muted"], markeredgewidth=0,
           label="stratified outlier, uncorroborated"),
    Line2D([0], [0], marker="o", linestyle="", markersize=5.5,
           markerfacecolor=C["warm"], markeredgewidth=0,
           label=f"Tier-2 corroborated window (n={len(candidates)})"),
    Line2D([0], [0], marker="D", linestyle="", markersize=4.5,
           markerfacecolor=HARD_Z, markeredgewidth=0,
           label=f"chrZ Tier-1, diploSHIC-only (n={len(zt1)})"),
]
ax.legend(handles=legend_handles, loc="lower left", bbox_to_anchor=(0.0, 1.01),
          ncol=4, handletextpad=0.3, columnspacing=1.1, borderaxespad=0, fontsize=6.5)

# x-axis: centered chromosome ticks
tick_pos = [offset[c] + chrom_len[c] / 2 for c in chroms]
tick_lab = [(c if (c.isdigit() and int(c) % 2 == 1) else ("Z" if c == "Z" else "")) for c in chroms]
ax.set_xticks(tick_pos)
ax.set_xticklabels(tick_lab, fontsize=6)
ax.set_xlim(-cum * 0.006, cum * 1.006)   # margins so chrZ markers are not clipped
ymax = max(r["y"] for r in rows)
ax.set_ylim(-0.15, ymax + 0.5)
ax.set_yticks([0, Y_THRESH, ymax])
ax.set_yticklabels(["0", "0.99", "≈1"])
ax.set_xlabel("chromosome")
ax.set_ylabel("stratified S percentile\n" r"$-\log_{10}(1-\mathrm{S\_pct})$")

despine(ax)
fig.tight_layout()
save(fig, OUT)
