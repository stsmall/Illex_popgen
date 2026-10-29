"""Supplementary Fig. S14. How the chromosome-2 inversion trajectory was fitted by grid search.

(a) the four fitted parameters as a frequency trajectory, with rejected one-parameter alternatives,
(b) the scoring recipe, (c, d) Delta-chi2 over the grid, (e) the four profiles with the parametric-bootstrap
cutoffs, (f) the per-statistic chi2 contributions that exclude each neighbouring grid value.

Reads the rebuilt grid (results/illex/refit_decline_pboot*_reps.jsonl.gz) and the bootstrap
(results/illex/param_bootstrap.json) from the msinv repository at commit 35412bf, and imports illex from it.
Run from an msinv checkout with its virtualenv:
  cd <msinv checkout> && PYTHONPATH=$PWD .venv/bin/python <manuscript>/scripts/plot_inversion_grid.py
"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figstyle import apply, C
apply()
from illex.balancing import decline_curve
from illex.scripts.param_bootstrap import load_cells, stats, profiles, TARGETS, SE, PARS

OUT = Path(__file__).resolve().parent.parent / "supp_figures" / "figS14_inversion_grid.png"
STAT_LAB = [r"$\pi_I/\pi_S$", r"$d_{xy}/\pi_I$", "singleton ratio"]
STAT_COL = ["#E69F00", "#56B4E9", "#009E73"]          # Okabe-Ito, validated
P_NOW = 0.374

cells = load_cells()
keys = np.array(sorted(cells))
means = np.array([stats(np.full(len(cells[tuple(k)][0]), 1 / len(cells[tuple(k)][0])), cells[tuple(k)])
                  for k in keys])
z2 = ((means - TARGETS) / SE) ** 2                 # per-statistic chi2 contributions
chi = z2.sum(1)
best = keys[chi.argmin()]
dchi = chi - chi.min()
prof = profiles(keys, chi)
boot = json.load(open("results/illex/param_bootstrap.json"))["params"]

fig = plt.figure(figsize=(13, 16))
gs = fig.add_gridspec(4, 4, height_ratios=[1.05, 1, 0.9, 1.0], hspace=0.6, wspace=0.42)

# ---------------- (a) what the four parameters mean ----------------
ax = fig.add_subplot(gs[0, :3])
# alternatives that change ONE parameter from the best fit; label placed where each curve departs
alts = [((850e3, 0.82, 175e3, 100e3), "plateau 0.82", (560, 0.835, "center")),
        ((850e3, 0.74, 0.0, 100e3), "fall ends today", (108, 0.64, "left")),
        ((850e3, 0.74, 100e3, 600e3), "600-ky fall", (455, 0.56, "center")),
        ((700e3, 0.74, 175e3, 100e3), "barrier 700 ky", (712, 0.80, "left"))]
for (ti, ph, td, tf), lab, (lx, ly, ha) in alts:
    t, p, _ = decline_curve(ti, td, tf, ph, P_NOW)
    k = np.flatnonzero((keys == (ti, ph, td, tf)).all(1))
    d = dchi[k[0]] if len(k) else np.nan
    ax.plot(np.asarray(t) / 1e3, p, color=C["muted"], lw=1.0, ls="--")
    ax.text(lx, ly, f"{lab}, " + r"$\Delta\chi^2$" + f" {d:.0f}", fontsize=7.5, color=C["muted"],
            va="center", ha=ha)
t, p, _ = decline_curve(best[0], best[2], best[3], best[1], P_NOW)
t = np.asarray(t) / 1e3
ax.plot(t, p, color=C["accent"], lw=2.4, label="best fit")
ti, ph, td, tf = best[0] / 1e3, best[1], best[2] / 1e3, best[3] / 1e3
y0 = 0.305
ax.annotate("", (ti, y0), (0, y0), arrowprops=dict(arrowstyle="<->", color=C["ink"], lw=0.9))
ax.text(ti / 2, y0 + 0.008, r"$t_\mathrm{inv}$ = age of the recombination barrier (850 ky)", ha="center",
        va="bottom", fontsize=8.5)
ax.annotate("", (td, 0.79), (td + tf, 0.79), arrowprops=dict(arrowstyle="<->", color=C["warm"], lw=0.9))
ax.text(td + tf / 2, 0.80, r"$t_\mathrm{fall}$ (100 ky)", ha="center", va="bottom", fontsize=8.5, color=C["warm"])
ax.annotate("", (0, 0.47), (td, 0.47), arrowprops=dict(arrowstyle="<->", color=C["warm"], lw=0.9))
ax.text(td / 2, 0.48, r"$t_\mathrm{decline}$ (175 ka)", ha="center", va="bottom", fontsize=8.5, color=C["warm"])
ax.text(560, ph - 0.015, r"$p_\mathrm{hist}$ = plateau (0.74), best fit", ha="center", va="top", fontsize=8.5,
        color=C["accent"])
ax.axhline(P_NOW, color=C["faint"], lw=1, zorder=0)
ax.text(995, P_NOW + 0.008, "today 0.374", fontsize=7.5, color=C["muted"], ha="right")
ax.set_xlim(0, 1000); ax.set_ylim(0.28, 0.87)
ax.set_xlabel("time before present (ky)  ←  present          past  →")
ax.set_ylabel("inverted-arrangement frequency")
ax.set_title("a  Each grid cell is a frequency trajectory (grey, alternatives that change one parameter)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---------------- (b) recipe ----------------
ax = fig.add_subplot(gs[0, 3]); ax.axis("off")
txt = ("How each grid cell is scored\n\n"
       "1. Pick one value per parameter\n    (one grid cell, 202 cells).\n"
       "2. Simulate the inversion 48 times\n    with msinv under that trajectory.\n"
       "3. Average three statistics over\n    the 48 replicates.\n"
       "4. Score the cell against the data\n"
       r"    $\chi^2=\sum_k\,[(\mathrm{sim}_k-\mathrm{obs}_k)/\mathrm{SE}_k]^2$" "\n"
       "5. Best cell = smallest " r"$\chi^2$" ".\n"
       "6. Profile: for each value of one\n    parameter, take the best cell over\n"
       "    the other three, report " r"$\Delta\chi^2$" ".\n"
       "7. A parametric bootstrap sets\n    the cutoff for exclusion (e).")
ax.text(0, 1, txt, va="top", fontsize=8.6, linespacing=1.35, color=C["ink"])
ax.set_title("b  Scoring a grid cell", loc="left")


# ---------------- (c, d) chi2 maps ----------------
def heat(ax, i, j, title):
    xs, ys = np.unique(keys[:, i]), np.unique(keys[:, j])
    M = np.full((len(ys), len(xs)), np.nan)
    for (xi, x) in enumerate(xs):
        for (yi, y) in enumerate(ys):
            m = (keys[:, i] == x) & (keys[:, j] == y)
            if m.any():
                M[yi, xi] = dchi[m].min()
    im = ax.imshow(M + 0.1, origin="lower", aspect="auto", cmap="Blues_r", norm=LogNorm(0.1, 120))
    for (xi, x) in enumerate(xs):
        for (yi, y) in enumerate(ys):
            v = M[yi, xi]
            ax.text(xi, yi, "not run" if np.isnan(v) else f"{v:.1f}", ha="center", va="center", fontsize=7.5,
                    color=C["muted"] if np.isnan(v) else ("white" if v < 3 else C["ink"]))
    return xs, ys, im


ax = fig.add_subplot(gs[1, :2])
xs, ys, im = heat(ax, 0, 1, "")
ax.set_xticks(range(len(xs)), [f"{x / 1e3:.0f}" for x in xs]); ax.set_yticks(range(len(ys)), [f"{y:.2f}" for y in ys])
ax.set_xlabel(r"$t_\mathrm{inv}$ (ky)"); ax.set_ylabel(r"$p_\mathrm{hist}$")
ax.set_title(r"c  $\Delta\chi^2$, age × plateau" "\n" "(fall timing re-optimised)")
ax = fig.add_subplot(gs[1, 2:])
xs, ys, im = heat(ax, 3, 2, "")
ax.set_xticks(range(len(xs)), [f"{x / 1e3:.0f}" for x in xs]); ax.set_yticks(range(len(ys)), [f"{y / 1e3:.0f}" for y in ys])
ax.set_xlabel(r"$t_\mathrm{fall}$ (ky)"); ax.set_ylabel(r"$t_\mathrm{decline}$ (ka)")
ax.set_title(r"d  $\Delta\chi^2$, fall duration × end of fall" "\n" "(age and plateau re-optimised)")
cb = fig.colorbar(im, ax=fig.axes[-2:], fraction=0.025, pad=0.02)
cb.set_label(r"$\Delta\chi^2$ from the best fit (darker is better)")

# ---------------- (e) profiles with calibrated cutoffs ----------------
lab = {"t_inv": (r"$t_\mathrm{inv}$ (ky)", 1e3), "p_hist": (r"$p_\mathrm{hist}$", 1),
       "t_decline": (r"$t_\mathrm{decline}$ (ka)", 1e3), "t_fall": (r"$t_\mathrm{fall}$ (ky)", 1e3)}
for j, p in enumerate(PARS):
    ax = fig.add_subplot(gs[2, j])
    v = np.array(sorted(prof[p])); d = np.array([prof[p][x] for x in v])
    c95 = boot[p]["c95"]
    ok = d <= c95
    ax.axhspan(0, c95, color=C["faint"], alpha=0.6, lw=0)
    ax.plot(v / lab[p][1], d + 0.1, color=C["ink"], lw=1.2)
    ax.scatter(v[ok] / lab[p][1], d[ok] + 0.1, s=30, color=C["accent"], zorder=3, label="inside 95% set")
    ax.scatter(v[~ok] / lab[p][1], d[~ok] + 0.1, s=30, facecolor="white", edgecolor=C["ink"], zorder=3,
               label="excluded")
    ax.axhline(2, color=C["muted"], ls=":", lw=1)
    ax.axhline(c95, color=C["warm"], ls="--", lw=1)
    ax.set_yscale("log"); ax.set_ylim(0.08, 150)
    ax.set_xlabel(lab[p][0])
    ax.text(0.02, 0.97, f"bootstrap cutoff {c95:.1f}", transform=ax.transAxes, fontsize=7.5, color=C["warm"],
            va="top")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if j == 2:
        ax.legend(loc="center right", fontsize=7)
    if j == 0:
        ax.set_ylabel(r"profile $\Delta\chi^2$ (+0.1, log)")
        ax.set_title(r"e  Profiles, shaded below the bootstrap-calibrated cutoff (dashed), dotted line at $\Delta\chi^2=2$",
                     loc="left")
        pass

# ---------------- (f) why the neighbours are rejected ----------------
ax = fig.add_subplot(gs[3, :])
def cell(ti, ph, td, tf):
    """best cell with these values for the named parameters (None = re-optimise)."""
    m = np.ones(len(keys), bool)
    for j, v in enumerate((ti, ph, td, tf)):
        if v is not None:
            m &= keys[:, j] == v
    return np.flatnonzero(m)[chi[m].argmin()]
rows = [("best fit (850 ky, 0.74, 175 ka, 100 ky)", cell(*best)),
        ("barrier 750 ky", cell(750e3, None, None, None)),
        ("barrier 900 ky", cell(900e3, None, None, None)),
        ("barrier 700 ky", cell(700e3, None, None, None)),
        ("plateau 0.78", cell(None, 0.78, None, None)),
        ("fall took 400 ky", cell(None, None, None, 400e3)),
        ("fall ended 50 ka", cell(None, None, 50e3, None)),
        ("fall ends today", cell(None, None, 0.0, None))]
y = np.arange(len(rows))[::-1]
for yi, (name, k) in zip(y, rows):
    left = 0.0
    for s_, col in enumerate(STAT_COL):
        w = z2[k, s_]
        ax.barh(yi, w, left=left, height=0.62, color=col, edgecolor="white", lw=1.5)
        if w > 1.2:
            ax.text(left + w / 2, yi, f"{w:.1f}", ha="center", va="center", fontsize=7.5, color=C["ink"])
        left += w
    ax.text(left + 0.6, yi, r"$\chi^2$ " + f"{left:.1f}", va="center", fontsize=8, color=C["ink"])
ax.set_yticks(y, [r[0] for r in rows])
ax.set_xlabel(r"contribution to $\chi^2$  $[(\mathrm{sim}-\mathrm{obs})/\mathrm{SE}]^2$, other parameters re-optimised")
ax.set_xlim(0, max(chi[r[1]] for r in rows) * 1.12)
ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(STAT_COL, STAT_LAB)], loc="upper right",
          ncol=3, fontsize=8)
ax.set_title("f  Which statistic rules out each neighbouring value", loc="left")
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)

fig.savefig(OUT, dpi=170, bbox_inches="tight")
print("wrote", OUT)
