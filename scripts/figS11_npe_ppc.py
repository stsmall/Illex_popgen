"""Supp Fig S11 -- posterior predictive checks of the balancing model (NPE draws).

For each statistic, the quantile of the observed value among 1,000 replicates simulated from NPE posterior draws of
the balancing model with gene flux, the favoured model, simulated (a) as 40-kb segments like the training
simulations and (b) as 474 independent 40-kb segments (18.96 Mb, the length of the inversion). Background bands mark the central 50, 90 and 95%
of the predictive distribution. Points are coloured by how far into the tails the observed value falls, and every
statistic with a two-sided predictive p below 0.10 is labelled with its p. Statistics withheld from the fit are
marked. The title gives the joint predictive p.

Data for (b) are ppc_per_stat.csv and ppc_joint.csv from illex/abc/npe/ppc_blocks.py + collect_ppc_blocks.py
(msinv results/illex/ppc_blocks_scored), copied to data/inversion_abc/npe_ppc_full_*.csv.
Data for (a) (written by illex/abc/smc_results/ppc/ppc_analyze.py on the output of illex/abc/npe/ppc_npe.py in the msinv
repository, branch illex/trajectory-abc): data/inversion_abc/npe_ppc_per_stat.csv, npe_ppc_joint.csv.
Output: supp_figures/figS11_npe_ppc.png
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from figstyle import apply, C, despine, save  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

apply()
D = ROOT / "data" / "inversion_abc"
SETS = [("a", "40-kb segments", pd.read_csv(D / "npe_ppc_per_stat.csv"), pd.read_csv(D / "npe_ppc_joint.csv")),
        ("b", "full length (474 × 40 kb)", pd.read_csv(D / "npe_ppc_full_per_stat.csv"),
         pd.read_csv(D / "npe_ppc_full_joint.csv"))]

LD_BINS = ["0.1–0.24", "0.24–0.58", "0.58–1.4", "1.4–3.4", "3.4–8.3", "8.3–20"]
GROUPS = [
    ("Diversity and divergence", [("pi_ratio", "π 2a / π 2+"), ("dxy_pi", "dxy / π 2a"),
                                  ("f1_ratio", "singleton ratio")]),
    ("Withheld from the fit", [("tajd_I", "Tajima's D, 2a"), ("tajd_S", "Tajima's D, 2+"),
                               ("edge_centre", "dxy/π edge / centre")]),
    ("Site-frequency spectrum, 2a", [(f"sfs_I_{b}", f"bin {b.replace('p', '+')}") for b in
                                     ("1", "2", "3", "4", "5", "6p")]),
    ("Site-frequency spectrum, 2+", [(f"sfs_S_{b}", f"bin {b.replace('p', '+')}") for b in
                                     ("1", "2", "3", "4", "5", "6p")]),
    ("LD (r²) by distance, 2a", [(f"ld_I_{k}", f"{LD_BINS[k - 1]} kb") for k in range(1, 7)]),
    ("LD (r²) by distance, 2+", [(f"ld_S_{k}", f"{LD_BINS[k - 1]} kb") for k in range(1, 7)]),
    ("Joint spectrum (2a bin, 2+ bin)", [(f"jsfs_{c}", c[0] + "," + c[1]) for c in
                                          ("01", "02", "03", "04", "10", "11", "12", "13", "14", "20", "21", "22")]),
]

# rows top to bottom, with a gap row before each group header
rows, headers = [], []
y = 0.0
for title, items in GROUPS:
    headers.append((y, title))
    y += 1.0
    for key, lab in items:
        rows.append((y, key, lab))
        y += 1.0
    y += 0.4
ymax = y

TAIL = [(0.05, "#D55E00", "$p<0.05$"), (0.10, "#E69F00", "$p$ 0.05–0.10"), (1.01, "#0072B2", "$p>0.10$")]


def colour(p):
    for lim, col, _ in TAIL:
        if p < lim:
            return col


FLUX = 1
fig, axs = plt.subplots(1, 2, figsize=(8.6, 11.2), sharey=True)
for ax, (letter, sub, per, joint) in zip(axs, SETS):
    d = per[per.flux == FLUX].set_index("stat")
    jp = float(joint.loc[joint.flux == FLUX, "joint_ppp"].iloc[0])
    for lo, hi, col in ((0, 0.025, "#bdbdbd"), (0.975, 1, "#bdbdbd"), (0.025, 0.05, "#dedede"),
                        (0.95, 0.975, "#dedede"), (0.05, 0.25, "#f2f2f2"), (0.75, 0.95, "#f2f2f2")):
        ax.axvspan(lo, hi, color=col, lw=0, zorder=0)
    ax.axvline(0.5, color=C["muted"], lw=0.7, ls=":", zorder=1)
    for yy, key, lab in rows:
        if key not in d.index:
            continue
        q, p = float(d.loc[key, "obs_quantile"]), float(d.loc[key, "ppp_two_sided"])
        ax.plot([0.5, q], [yy, yy], color=C["faint"], lw=1.0, zorder=2)
        ax.scatter(q, yy, s=34, color=colour(p), edgecolor="white", lw=0.5, zorder=3)
        if p < 0.10:
            ha, dx = ("right", -0.025) if q > 0.5 else ("left", 0.025)
            ax.text(q + dx, yy, "p < 0.01" if p < 0.005 else f"p = {p:.2f}", ha=ha, va="center", fontsize=7.5, color=C["ink"], zorder=4)
    for yy, title in headers:
        ax.text(0.0, yy + 0.15, title, fontsize=8.5, fontweight="bold", color=C["ink"], va="center", ha="left")
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.025, 0.25, 0.5, 0.75, 0.975, 1])
    ax.set_xticklabels(["0", "", "0.25", "0.5", "0.75", "", "1"])
    ax.set_xlabel("Quantile of observed value\namong predictive replicates")
    ax.set_title(f"{sub}\njoint predictive p = {jp:.2f}", fontsize=10)
    ax.text(-0.02, 1.035, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", ha="right")
    despine(ax, keep=("bottom",))
    ax.tick_params(axis="y", length=0)
axs[0].set_yticks([r[0] for r in rows])
axs[0].set_yticklabels([r[2] for r in rows], fontsize=8)
axs[0].set_ylim(ymax - 0.2, -0.6)
handles = [plt.Line2D([], [], marker="o", ls="", color=col, markersize=6, label=lab) for _, col, lab in TAIL[::-1]]
handles += [plt.Rectangle((0, 0), 1, 1, color=c, label=l) for c, l in
            (("#f2f2f2", "central 50–90%"), ("#dedede", "90–95%"), ("#bdbdbd", "outer 5%"))]
fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.6, -0.01),
           handletextpad=0.4, columnspacing=1.0)
fig.subplots_adjust(left=0.2, right=0.98, top=0.94, bottom=0.09, wspace=0.08)
save(fig, str(ROOT / "supp_figures" / "figS11_npe_ppc.png"))
