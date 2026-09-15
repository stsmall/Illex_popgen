"""DRAFT recombination-summary ideas (recomb as a summary of the evolutionary process
and its link to Ne / selection). Writes several candidate figures to
manuscript/figures/recomb_draft_figs/ for the user to choose among.

Data: steps/14_sweep_seqmodel/results/bgs_diagnostic_windows.tsv
      (per ~50 kb window: pi, tajD, cds_frac, recomb [M/bp]) + /tmp/win_cds_dist.tsv
      (window -> distance to nearest CDS, bp, by row order).
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats

OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/recomb_draft_figs"
os.makedirs(OUT, exist_ok=True)
W = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/bgs_diagnostic_windows.tsv"

d = pd.read_csv(W, sep="\t")
d["chrom"] = d["chrom"].astype(str)
d["cMMb"] = d["recomb"].astype(float) * 1e8            # M/bp -> cM/Mb
dist = pd.read_csv("/tmp/win_cds_dist.tsv", sep="\t", header=None, names=["row", "dist_cds"])
dist = dist.groupby("row", as_index=True)["dist_cds"].min()   # dedupe bedtools ties
d = d.reset_index(drop=True)
d["dist_cds"] = dist.reindex(np.arange(1, len(d) + 1)).to_numpy()
d = d[np.isfinite(d["cMMb"]) & np.isfinite(d["pi"])].copy()
ACC = "#0072B2"

def binmean(x, y, q=12):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    edges = np.quantile(x, np.linspace(0, 1, q + 1))
    edges = np.unique(edges)
    idx = np.clip(np.digitize(x, edges[1:-1]), 0, len(edges) - 2)
    mids = 0.5 * (edges[:-1] + edges[1:])
    m = np.array([y[idx == i].mean() for i in range(len(mids))])
    se = np.array([y[idx == i].std() / max(np.sqrt((idx == i).sum()), 1) for i in range(len(mids))])
    return mids, m, se

# ---------- Draft 1: recombination landscape vs relative chromosome position ----------
d["relpos"] = d.groupby("chrom")["start"].transform(lambda s: (s - s.min()) / (s.max() - s.min() + 1))
fig, ax = plt.subplots(figsize=(6.6, 4.0))
mids, m, se = binmean(d["relpos"].to_numpy(), d["cMMb"].to_numpy(), q=25)
ax.fill_between(mids, m - se, m + se, color=ACC, alpha=0.2, lw=0)
ax.plot(mids, m, color=ACC, lw=2)
ax.set_xlabel("relative position along chromosome (0 = start, 1 = end)")
ax.set_ylabel("recombination rate (cM/Mb)")
ax.set_title("Draft 1 — recombination landscape (mean across autosomes)", loc="left", fontweight="bold")
despine(ax); fig.savefig(f"{OUT}/draft1_landscape.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------- Draft 2: recombination vs genic architecture ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={"wspace": 0.3})
# (a) recomb vs distance to nearest CDS
dd = d[np.isfinite(d["dist_cds"])]
mids, m, se = binmean(np.log10(dd["dist_cds"].to_numpy() + 1), dd["cMMb"].to_numpy(), q=14)
ax[0].fill_between(mids, m - se, m + se, color=ACC, alpha=0.2, lw=0)
ax[0].plot(mids, m, color=ACC, lw=2, marker="o", ms=3)
ax[0].set_xlabel(r"$\log_{10}$ distance to nearest exon/CDS (bp)")
ax[0].set_ylabel("recombination rate (cM/Mb)")
ax[0].set_title("a  recomb vs distance from genes", loc="left")
despine(ax[0])
# (b) genic vs intergenic + CDS-density bins
genic = d[d["cds_frac"] > 0]["cMMb"]; inter = d[d["cds_frac"] == 0]["cMMb"]
mids, m, se = binmean(d[d["cds_frac"] > 0]["cds_frac"].to_numpy(), d[d["cds_frac"] > 0]["cMMb"].to_numpy(), q=12)
ax[1].fill_between(mids, m - se, m + se, color="#009E73", alpha=0.2, lw=0)
ax[1].plot(mids, m, color="#009E73", lw=2, marker="o", ms=3)
ax[1].axhline(inter.mean(), ls="--", color=C["muted"], lw=1)
ax[1].text(mids.max(), inter.mean(), "  intergenic mean", fontsize=7.5, color=C["muted"], va="center")
ax[1].set_xlabel("local CDS density (fraction coding)")
ax[1].set_ylabel("recombination rate (cM/Mb)")
ax[1].set_title("b  recomb vs coding density", loc="left")
despine(ax[1])
fig.suptitle("Draft 2 — recombination and genic architecture", x=0.02, ha="left", fontweight="bold")
fig.savefig(f"{OUT}/draft2_genic.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------- Draft 3: the selection / Ne connection (Begun-Aquadro) ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={"wspace": 0.3})
r = d["cMMb"].to_numpy()
# (a) pi vs recombination
mids, m, se = binmean(r, d["pi"].to_numpy(), q=14)
rho = stats.spearmanr(r, d["pi"], nan_policy="omit").correlation
ax[0].fill_between(mids, m - se, m + se, color=ACC, alpha=0.2, lw=0)
ax[0].plot(mids, m, color=ACC, lw=2, marker="o", ms=3)
ax[0].set_xlabel("recombination rate (cM/Mb)")
ax[0].set_ylabel(r"nucleotide diversity $\pi$")
ax[0].set_title(f"a  diversity vs recombination (Spearman $\\rho$={rho:+.2f})", loc="left")
ax[0].text(0.5, 0.06, "positive slope = linked selection (BGS + sweeps)\nreduces diversity where recombination is low\n$\\Rightarrow$ locally reduced effective $N_e$",
           transform=ax[0].transAxes, fontsize=7.5, color=C["muted"], va="bottom")
despine(ax[0])
# (b) Tajima's D vs recombination
mids, m, se = binmean(r, d["tajD"].to_numpy(), q=14)
rho2 = stats.spearmanr(r, d["tajD"], nan_policy="omit").correlation
ax[1].fill_between(mids, m - se, m + se, color="#D55E00", alpha=0.2, lw=0)
ax[1].plot(mids, m, color="#D55E00", lw=2, marker="o", ms=3)
ax[1].set_xlabel("recombination rate (cM/Mb)")
ax[1].set_ylabel(r"Tajima's $D$")
ax[1].set_title(f"b  Tajima's $D$ vs recombination (Spearman $\\rho$={rho2:+.2f})", loc="left")
despine(ax[1])
fig.suptitle("Draft 3 — recombination, diversity and the efficacy of selection", x=0.02, ha="left", fontweight="bold")
fig.savefig(f"{OUT}/draft3_selection_Ne.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------- Draft 4: combined candidate summary (2x2) ----------
fig, axs = plt.subplots(2, 2, figsize=(10.5, 8.0), gridspec_kw={"wspace": 0.32, "hspace": 0.42})
mids, m, se = binmean(d["relpos"].to_numpy(), d["cMMb"].to_numpy(), q=25)
axs[0, 0].fill_between(mids, m - se, m + se, color=ACC, alpha=0.2, lw=0); axs[0, 0].plot(mids, m, color=ACC, lw=2)
axs[0, 0].set_xlabel("relative chromosome position (0–1)"); axs[0, 0].set_ylabel("cM/Mb")
axs[0, 0].set_title("a  recombination landscape", loc="left", fontweight="bold"); despine(axs[0, 0])
mids, m, se = binmean(np.log10(dd["dist_cds"].to_numpy() + 1), dd["cMMb"].to_numpy(), q=14)
axs[0, 1].fill_between(mids, m - se, m + se, color="#009E73", alpha=0.2, lw=0); axs[0, 1].plot(mids, m, color="#009E73", lw=2, marker="o", ms=3)
axs[0, 1].set_xlabel(r"$\log_{10}$ dist. to exon (bp)"); axs[0, 1].set_ylabel("cM/Mb")
axs[0, 1].set_title("b  recomb vs distance from genes", loc="left", fontweight="bold"); despine(axs[0, 1])
mids, m, se = binmean(r, d["pi"].to_numpy(), q=14)
axs[1, 0].fill_between(mids, m - se, m + se, color=ACC, alpha=0.2, lw=0); axs[1, 0].plot(mids, m, color=ACC, lw=2, marker="o", ms=3)
axs[1, 0].set_xlabel("recombination rate (cM/Mb)"); axs[1, 0].set_ylabel(r"$\pi$")
axs[1, 0].set_title(f"c  diversity vs recomb ($\\rho$={rho:+.2f})", loc="left", fontweight="bold"); despine(axs[1, 0])
mids, m, se = binmean(r, d["tajD"].to_numpy(), q=14)
axs[1, 1].fill_between(mids, m - se, m + se, color="#D55E00", alpha=0.2, lw=0); axs[1, 1].plot(mids, m, color="#D55E00", lw=2, marker="o", ms=3)
axs[1, 1].set_xlabel("recombination rate (cM/Mb)"); axs[1, 1].set_ylabel(r"Tajima's $D$")
axs[1, 1].set_title(f"d  Tajima's $D$ vs recomb ($\\rho$={rho2:+.2f})", loc="left", fontweight="bold"); despine(axs[1, 1])
fig.suptitle("Draft 4 — combined: recombination as a summary of linked selection", x=0.02, ha="left", fontweight="bold", fontsize=13)
fig.savefig(f"{OUT}/draft4_combined.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------- Draft 5: RECOMMENDED -- Ne/selection via the INDEPENDENT gene-density axis ----------
# (avoids the ReLERNN circularity: recomb is inferred FROM diversity/LD, so pi~recomb is
#  not an independent test. Gene/CDS density is independent, and is where illex's linked
#  selection actually lives -- BGS is gene-density-driven, not recombination-driven.)
rho_cds = stats.spearmanr(d["cds_frac"], d["pi"], nan_policy="omit").correlation
ldist = np.log10(dd["dist_cds"].to_numpy() + 1)
fig, axs = plt.subplots(2, 2, figsize=(11.0, 8.0), gridspec_kw={"wspace": 0.30, "hspace": 0.40})
# (a) recomb vs distance from exon
mids, m, se = binmean(ldist, dd["cMMb"].to_numpy(), q=14)
axs[0, 0].fill_between(mids, m - se, m + se, color="#009E73", alpha=0.2, lw=0)
axs[0, 0].plot(mids, m, color="#009E73", lw=2, marker="o", ms=3)
axs[0, 0].set_xlabel(r"$\log_{10}$ distance to nearest exon (bp)"); axs[0, 0].set_ylabel("recombination (cM/Mb)")
axs[0, 0].set_title("a  recombination rises away from genes", loc="left", fontweight="bold"); despine(axs[0, 0])
# (b) pi vs distance from exon (MATCHING panel a)
mids, m, se = binmean(ldist, dd["pi"].to_numpy(), q=14)
rho_pd = stats.spearmanr(ldist, dd["pi"], nan_policy="omit").correlation
axs[0, 1].fill_between(mids, m - se, m + se, color="#0072B2", alpha=0.2, lw=0)
axs[0, 1].plot(mids, m, color="#0072B2", lw=2, marker="o", ms=3)
axs[0, 1].set_xlabel(r"$\log_{10}$ distance to nearest exon (bp)"); axs[0, 1].set_ylabel(r"nucleotide diversity $\pi$")
axs[0, 1].set_title(f"b  diversity also rises away from genes ($\\rho$={rho_pd:+.2f})", loc="left", fontweight="bold")
axs[0, 1].text(0.5, 0.08, "both recover with distance from coding\nsequence — the footprint of linked selection",
               transform=axs[0, 1].transAxes, fontsize=7.5, color=C["muted"], va="bottom")
despine(axs[0, 1])
# (c) pi vs CDS density -- the Ne/linked-selection signal (independent of recomb)
mids, m, se = binmean(d["cds_frac"].to_numpy(), d["pi"].to_numpy(), q=14)
axs[1, 0].fill_between(mids, m - se, m + se, color="#0072B2", alpha=0.2, lw=0)
axs[1, 0].plot(mids, m, color="#0072B2", lw=2, marker="o", ms=3)
axs[1, 0].set_xlabel("local coding density (fraction CDS)"); axs[1, 0].set_ylabel(r"nucleotide diversity $\pi$")
axs[1, 0].set_title(f"c  diversity falls in gene-dense regions ($\\rho$={rho_cds:+.2f})", loc="left", fontweight="bold")
axs[1, 0].text(0.5, 0.9, "background selection reduces $\\pi$ where\ngenes are dense $\\Rightarrow$ locally reduced $N_e$",
               transform=axs[1, 0].transAxes, fontsize=7.5, color=C["muted"], va="top")
despine(axs[1, 0])
# (d) pi vs recomb -- shown as the ABSENCE of the classic relationship
mids, m, se = binmean(r, d["pi"].to_numpy(), q=14)
axs[1, 1].fill_between(mids, m - se, m + se, color=C["muted"], alpha=0.2, lw=0)
axs[1, 1].plot(mids, m, color=C["muted"], lw=2, marker="o", ms=3)
axs[1, 1].set_xlabel("recombination rate (cM/Mb)"); axs[1, 1].set_ylabel(r"nucleotide diversity $\pi$")
axs[1, 1].set_title(f"d  no diversity–recombination relation ($\\rho$={rho:+.2f})", loc="left", fontweight="bold")
axs[1, 1].text(0.5, 0.9, "unlike Drosophila/humans (and ReLERNN\ninfers recomb FROM diversity — not independent)",
               transform=axs[1, 1].transAxes, fontsize=7.0, color=C["muted"], va="top")
despine(axs[1, 1])
fig.suptitle("Draft 5 (recommended) — gene proximity shapes both recombination and diversity",
             x=0.02, ha="left", fontweight="bold", fontsize=13)
fig.savefig(f"{OUT}/draft5_recommended.png", dpi=200, bbox_inches="tight"); plt.close(fig)

print("wrote 5 drafts to", OUT)
print(f"Spearman pi~recomb = {rho:+.3f}; tajD~recomb = {rho2:+.3f}; pi~cds_density = {rho_cds:+.3f}; n windows = {len(d)}")
