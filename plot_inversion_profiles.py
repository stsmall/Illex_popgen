"""Supp: the two profile-likelihood scans that bound the chr2 inversion's decline
(moved out of Fig 3). b-ii: when the fall ended (t_decline); b-iii: how long it took (t_fall).
Delta-chi-square from the best fit; the shaded band is Delta-chi2 < 2 (not excluded).
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/dev/shm/mplcache")
import csv, sys
from pathlib import Path
sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, OKABE, despine, save
apply()
import numpy as np, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

RESULTS_ILLEX = Path("/sietch_colab/ssmall/projects/msinv_dir/inversion_sims/files/"
                     "results/illex")
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/supp_figures/supp_inversion_profiles.png"
COL_BB = OKABE[7]
COL_GONE2 = OKABE[3]

TARGETS = {"pi_i_over_pi_s": (1.3556, 0.0481),
           "dxy_over_pi_i":  (1.3848, 0.0214),
           "f1_ratio":       (0.8125, 0.0076)}


def chi2(row):
    return sum(((float(row[k]) - t) / se) ** 2 for k, (t, se) in TARGETS.items())


def profile(tag, key):
    rows = list(csv.DictReader(open(RESULTS_ILLEX / f"refit_decline{tag}.csv")))
    best = {}
    for r in rows:
        v, c = float(r[key]), chi2(r)
        if v not in best or c < best[v][0]:
            best[v] = (c, r)
    gmin = min(c for c, _ in best.values())
    vals = sorted(best)
    return np.array(vals), np.array([best[v][0] - gmin for v in vals])


fig = plt.figure(figsize=(7.6, 3.4))
gs = GridSpec(1, 2, figure=fig, wspace=0.28)
for j, (tag, key, sub, xlab, subtitle) in enumerate((
        ("tdecline", "t_decline", "a", r"$t_{decline}$  (ka)", "When the fall ended"),
        ("tfall2", "t_fall", "b", r"$t_{fall}$  (ky)", "How long the fall took"))):
    ax = fig.add_subplot(gs[0, j])
    v, dc = profile(tag, key)
    ax.axhspan(0, 2, color=COL_GONE2, alpha=0.10, lw=0)
    ax.axhline(2, lw=0.7, ls=":", color=C["muted"])
    ax.plot(v / 1e3, dc, marker="o", ms=4.2, lw=1.2, color=COL_BB,
            mfc="white", mec=COL_BB, mew=1.2, clip_on=False, zorder=3)
    ax.set_yscale("symlog", linthresh=2)
    ax.set_ylim(-0.4, 260)
    ax.set_yticks([0, 1, 2, 10, 100]); ax.set_yticklabels(["0", "1", "2", "10", "100"])
    ax.set_xlabel(xlab)
    if j == 0:
        ax.set_ylabel(r"$\Delta\chi^2$ from best fit")
    despine(ax)
    for x, y in zip(v / 1e3, dc):
        if y > 2:
            ax.annotate(f"{y:.0f}" if y >= 10 else f"{y:.1f}", xy=(x, y),
                        xytext=(0, 5), textcoords="offset points",
                        ha="center", fontsize=7.5, color=C["ink"])
    ax.text(-0.18, 1.08, sub, transform=ax.transAxes, fontweight="bold", fontsize=12)
    ax.set_title(subtitle, loc="left", pad=8, fontsize=9.5, fontweight="normal", color=C["muted"])

save(fig, OUT)
print("wrote", OUT)
