"""Supplement: recombination, gene proximity and the efficacy of selection in Illex.
Gene proximity shapes BOTH recombination and diversity (the footprint of linked selection);
diversity tracks gene density (background selection -> locally reduced Ne); and there is NO
diversity-recombination relationship (unlike Drosophila/humans; ReLERNN also infers recomb
from diversity, so the two are not independent).

Data (per ~50 kb window): steps/14_sweep_seqmodel/results/bgs_diagnostic_windows.tsv
      (pi, tajD, cds_frac, recomb [M/bp]) + .../win_cds_dist.tsv (dist to nearest CDS, bp).
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import sys
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/scripts")
from figstyle import apply, C, despine
apply()
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats

R = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/figS17_recomb_selection.png"

d = pd.read_csv(f"{R}/bgs_diagnostic_windows.tsv", sep="\t")
d["chrom"] = d["chrom"].astype(str)
d["cMMb"] = d["recomb"].astype(float) * 1e8
dist = pd.read_csv(f"{R}/win_cds_dist.tsv", sep="\t", header=None, names=["row", "dist_cds"])
dist = dist.groupby("row")["dist_cds"].min()
d = d.reset_index(drop=True)
d["dist_cds"] = dist.reindex(np.arange(1, len(d) + 1)).to_numpy()
d = d[np.isfinite(d["cMMb"]) & np.isfinite(d["pi"])].copy()
dd = d[np.isfinite(d["dist_cds"])]
ldist = np.log10(dd["dist_cds"].to_numpy() + 1)
r = d["cMMb"].to_numpy()
GREEN, BLUE, GREY = "#009E73", "#0072B2", C["muted"]


def binmean(x, y, q=14):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    edges = np.unique(np.quantile(x, np.linspace(0, 1, q + 1)))
    idx = np.clip(np.digitize(x, edges[1:-1]), 0, len(edges) - 2)
    mids = 0.5 * (edges[:-1] + edges[1:])
    m = np.array([y[idx == i].mean() for i in range(len(mids))])
    se = np.array([y[idx == i].std() / max(np.sqrt((idx == i).sum()), 1) for i in range(len(mids))])
    return mids, m, se


fig, axs = plt.subplots(2, 2, figsize=(10.8, 8.0), gridspec_kw={"wspace": 0.30, "hspace": 0.42})

# (a) recomb vs distance from exon -- the CONTROL. Drawn on the full recombination-rate
#     axis (same scale as panel d) so the near-flatness is visible: the whole gradient is
#     0.206 -> 0.217 cM/Mb (~6%), Spearman rho ~ +0.05, and it is unchanged when pi is held
#     fixed (partial rho +0.06) -- i.e. NOT a diversity/Ne leak; recombination is agnostic
#     to gene content, so it cannot be what organises the diversity gradients in (b)/(c).
mids, m, se = binmean(ldist, dd["cMMb"].to_numpy())
rho_rd = stats.spearmanr(ldist, dd["cMMb"], nan_policy="omit").correlation
near = dd.loc[dd["dist_cds"] < 1e3, "cMMb"].mean(); far = dd.loc[dd["dist_cds"] > 1e5, "cMMb"].mean()
axs[0, 0].fill_between(mids, m - se, m + se, color=GREEN, alpha=0.2, lw=0)
axs[0, 0].plot(mids, m, color=GREEN, lw=2, marker="o", ms=3)
axs[0, 0].set_ylim(0, max(0.32, r[np.isfinite(r)].mean() * 1.5))
axs[0, 0].set_xlabel(r"$\log_{10}$ distance to nearest exon (bp)")
axs[0, 0].set_ylabel("recombination rate (cM/Mb)")
axs[0, 0].set_title(f"a   recombination is flat vs gene proximity ($\\rho={rho_rd:+.2f}$)",
                    loc="left", fontweight="bold", fontsize=10)
axs[0, 0].text(0.04, 0.96,
               f"{near:.3f} → {far:.3f} cM/Mb near→far ({100*(far/near-1):.0f}% range);\n"
               "unchanged when $\\pi$ is held fixed — recombination\nis agnostic to gene content (the control)",
               transform=axs[0, 0].transAxes, fontsize=7.6, color=GREY, va="top", ha="left")
despine(axs[0, 0])

# (b) pi vs distance from exon (same x as a); annotation top-left (empty region)
mids, m, se = binmean(ldist, dd["pi"].to_numpy())
rho_pd = stats.spearmanr(ldist, dd["pi"], nan_policy="omit").correlation
axs[0, 1].fill_between(mids, m - se, m + se, color=BLUE, alpha=0.2, lw=0)
axs[0, 1].plot(mids, m, color=BLUE, lw=2, marker="o", ms=3)
axs[0, 1].set_xlabel(r"$\log_{10}$ distance to nearest exon (bp)")
axs[0, 1].set_ylabel(r"nucleotide diversity $\pi$")
axs[0, 1].set_title(f"b   $\\pi$ rises away from genes ($\\rho={rho_pd:+.2f}$)",
                    loc="left", fontweight="bold", fontsize=10)
axs[0, 1].text(0.04, 0.96, "$\\pi$ recovers with distance from coding\nsequence while recombination (a) does not —\nthe footprint of linked selection",
               transform=axs[0, 1].transAxes, fontsize=7.6, color=GREY, va="top", ha="left")
despine(axs[0, 1])

# (c) pi vs CDS density -- BGS / Ne signal (independent of recomb)
rho_cds = stats.spearmanr(d["cds_frac"], d["pi"], nan_policy="omit").correlation
mids, m, se = binmean(d["cds_frac"].to_numpy(), d["pi"].to_numpy())
axs[1, 0].fill_between(mids, m - se, m + se, color=BLUE, alpha=0.2, lw=0)
axs[1, 0].plot(mids, m, color=BLUE, lw=2, marker="o", ms=3)
axs[1, 0].set_xlabel("local coding density (fraction CDS)")
axs[1, 0].set_ylabel(r"nucleotide diversity $\pi$")
axs[1, 0].set_title(f"c   diversity falls in gene-dense regions ($\\rho={rho_cds:+.2f}$)",
                    loc="left", fontweight="bold")
axs[1, 0].text(0.5, 0.94, "background selection reduces $\\pi$ where\ngenes are dense $\\Rightarrow$ locally reduced $N_e$",
               transform=axs[1, 0].transAxes, fontsize=7.6, color=GREY, va="top", ha="left")
despine(axs[1, 0])

# (d) pi vs recomb -- the ABSENCE of the classic relationship
rho_r = stats.spearmanr(r, d["pi"], nan_policy="omit").correlation
mids, m, se = binmean(r, d["pi"].to_numpy())
axs[1, 1].fill_between(mids, m - se, m + se, color=GREY, alpha=0.25, lw=0)
axs[1, 1].plot(mids, m, color=GREY, lw=2, marker="o", ms=3)
axs[1, 1].set_xlabel("recombination rate (cM/Mb)")
axs[1, 1].set_ylabel(r"nucleotide diversity $\pi$")
axs[1, 1].set_title(f"d   no diversity–recombination relation ($\\rho={rho_r:+.2f}$)",
                    loc="left", fontweight="bold")
axs[1, 1].text(0.5, 0.96, "unlike Drosophila/humans: recombination is\nnearly uniform here, so it cannot organise diversity",
               transform=axs[1, 1].transAxes, fontsize=7.0, color=GREY, va="top", ha="left")
despine(axs[1, 1])

fig.savefig(OUT, dpi=220, bbox_inches="tight")
print("wrote", OUT)
print(f"rho: recomb~dist={rho_rd:+.3f} (near {near:.4f} far {far:.4f})  pi~dist={rho_pd:+.3f}  "
      f"pi~cds={rho_cds:+.3f}  pi~recomb={rho_r:+.3f}  n={len(d)}")
