"""Supplementary Fig. S14. How the chromosome-2 inversion trajectory was fitted by grid search.

(a) the four fitted parameters as a frequency trajectory, with the best trajectory at excluded values,
(b) the scoring recipe, (c, d) Delta D over the grid, (e) the four profiles with the parametric-bootstrap
cutoffs, (f) the sequential per-statistic contributions to D that exclude each neighbouring grid value.

The fit was run with the BB arrangement as the trajectory class (p_BB, a plateau then a fall to 0.374). The model is
symmetric in the arrangement labels, so everything is plotted for the derived AA arrangement, p_AA = 1 - p_BB
(a plateau then a rise to 0.626), with I = AA and S = BB in the statistic labels.

D is the covariance-aware distance (m - x)' S^-1 (m - x), with S the measurement covariance plus the trajectory
variance from the 825-ky best cell. Reads the grid (results/illex/refit_decline_pboot{TAGS}_reps.jsonl.gz) and
the recovery output (results/illex/param_recovery_{TAGS}_traj825.json) from the msinv repository
(branch illex/pboot-recovery), and imports illex from it.
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
STAT_LAB_F = [r"$\pi_S/\pi_I$", r"$d_{xy}/\pi_S$ given $\pi_S/\pi_I$", "singleton ratio given both"]
STAT_COL = ["#E69F00", "#56B4E9", "#009E73"]          # Okabe-Ito, validated
P_NOW = 0.374          # p_BB today; the plot shows p_AA = 1 - p_BB = 0.626

TAGS = sys.argv[1] if len(sys.argv) > 1 else "ABCDEFGHIJK"
cells = load_cells(TAGS)
keys = np.array(sorted(cells))
means = np.array([stats(np.full(len(cells[tuple(k)][0]), 1 / len(cells[tuple(k)][0])), cells[tuple(k)])
                  for k in keys])
rec = json.load(open(f"results/illex/param_recovery_{TAGS}_traj825.json"))
S = np.array(rec["total_cov"])
# sequential contributions: whiten the residuals with the Cholesky factor of S, so D = sum of squares and the
# k-th term is statistic k given the statistics before it
z2 = np.linalg.solve(np.linalg.cholesky(S), (means - TARGETS).T).T ** 2
chi = z2.sum(1)
best = keys[chi.argmin()]
dchi = chi - chi.min()
prof = profiles(keys, chi)
boot = rec["cov"]["params"]
assert np.allclose(best, [rec["cov"]["best"][p] for p in PARS])


def cell(ti, ph, td, tf):
    """best cell with these values for the named parameters (None = re-optimise)."""
    m = np.ones(len(keys), bool)
    for j, v in enumerate((ti, ph, td, tf)):
        if v is not None:
            m &= keys[:, j] == v
    return np.flatnonzero(m)[chi[m].argmin()]


def fmt(k):
    t, p, d, f = keys[k]
    return f"{t / 1e3:.0f} ky, {1 - p:.2f}, {d / 1e3:.0f} ka, {f / 1e3:.0f} ky"

fig = plt.figure(figsize=(13, 16))
gs = fig.add_gridspec(4, 4, height_ratios=[1.05, 1, 0.9, 1.0], hspace=0.6, wspace=0.42)

# ---------------- (a) what the four parameters mean ----------------
ax = fig.add_subplot(gs[0, :3])
# the best trajectory at an excluded value of one parameter (the other three re-optimised)
alts = [(cell(None, 0.82, None, None), "plateau 0.18", (560, 0.165, "center")),
        (cell(None, 0.56, None, None), "plateau 0.44", (560, 0.455, "center")),
        (cell(None, None, 0.0, None), "rise ends today", (8, 0.355, "left")),
        (cell(None, None, None, 600e3), "600-ky rise", (430, 0.55, "right")),
        (cell(700e3, None, None, None), "barrier 700 ky", (705, 0.235, "right"))]
for k, lab, (lx, ly, ha) in alts:
    ti, ph, td, tf = keys[k]
    t, p, _ = decline_curve(ti, td, tf, ph, P_NOW)
    ax.plot(np.asarray(t) / 1e3, 1 - np.asarray(p), color=C["muted"], lw=1.0, ls="--")
    ax.text(lx, ly, f"{lab}, " + r"$\Delta D$" + f" {dchi[k]:.0f}", fontsize=7.5, color=C["muted"],
            va="center", ha=ha)
t, p, _ = decline_curve(best[0], best[2], best[3], best[1], P_NOW)
t = np.asarray(t) / 1e3
ax.plot(t, 1 - np.asarray(p), color=C["accent"], lw=2.4, label="best fit")
ti, ph, td, tf = best[0] / 1e3, best[1], best[2] / 1e3, best[3] / 1e3
y0 = 0.665
ax.annotate("", (ti, y0), (0, y0), arrowprops=dict(arrowstyle="<->", color=C["ink"], lw=0.9))
ax.text(ti / 2, y0 + 0.008, r"$t_\mathrm{inv}$ = age of the recombination barrier" f" ({ti:.0f} ky)", ha="center",
        va="bottom", fontsize=8.5)
ax.annotate("", (td, 0.24), (td + tf, 0.24), arrowprops=dict(arrowstyle="<->", color=C["warm"], lw=0.9))
ax.text(td + tf / 2, 0.235, r"$t_\mathrm{rise}$" f" ({tf:.0f} ky)", ha="center", va="top", fontsize=8.5,
        color=C["warm"])
ax.annotate("", (0, 0.53), (td, 0.53), arrowprops=dict(arrowstyle="<->", color=C["warm"], lw=0.9))
ax.text(td / 2, 0.52, r"$t_\mathrm{end}$" f" ({td:.0f} ka)", ha="center", va="top", fontsize=8.5,
        color=C["warm"])
ax.text(560, 1 - ph + 0.06, r"$p_\mathrm{hist}$" f" = plateau ({1 - ph:.2f}), best fit", ha="center", va="bottom", fontsize=8.5,
        color=C["accent"])
ax.axhline(1 - P_NOW, color=C["faint"], lw=1, zorder=0)
ax.text(995, 1 - P_NOW + 0.008, "today 0.626", fontsize=7.5, color=C["muted"], ha="right")
ax.set_xlim(0, 1000); ax.set_ylim(0.12, 0.72)
ax.set_xlabel("time before present (ky)  ←  present          past  →")
ax.set_ylabel("derived (2a) arrangement frequency")
ax.set_title("a  Each grid cell is a frequency trajectory (grey, best trajectory at an excluded value)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---------------- (b) recipe ----------------
ax = fig.add_subplot(gs[0, 3]); ax.axis("off")
txt = ("How each grid cell is scored\n\n"
       f"1. Pick one value per parameter\n    (one grid cell, {len(keys):,} cells).\n"
       "2. Simulate the inversion 48 times\n    with msinv under that trajectory.\n"
       "3. Average three statistics over\n    the 48 replicates.\n"
       "4. Score the cell against the data\n"
       r"    $D=(\mathrm{sim}-\mathrm{obs})^\top\Sigma^{-1}(\mathrm{sim}-\mathrm{obs})$" "\n"
       "    " r"$\Sigma$" " = measurement covariance\n    + trajectory variance.\n"
       "5. Best cell = smallest " r"$D$" ".\n"
       "6. Profile: for each value of one\n    parameter, take the best cell over\n"
       "    the other three, report " r"$\Delta D$" ".\n"
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
    im = ax.imshow(M + 0.1, origin="lower", aspect="auto", cmap="Blues_r", norm=LogNorm(0.1, 150))
    for (xi, x) in enumerate(xs):
        for (yi, y) in enumerate(ys):
            v = M[yi, xi]
            ax.text(xi, yi, "–" if np.isnan(v) else f"{v:.1f}", ha="center", va="center", fontsize=6,
                    color=C["muted"] if np.isnan(v) else ("white" if v < 3 else C["ink"]))
    return xs, ys, im


ax = fig.add_subplot(gs[1, :2])
xs, ys, im = heat(ax, 0, 1, "")
ax.set_xticks(range(len(xs)), [f"{x / 1e3:.0f}" for x in xs]); ax.set_yticks(range(len(ys)), [f"{1 - y:.2f}" for y in ys])
ax.set_xlabel(r"$t_\mathrm{inv}$ (ky)"); ax.set_ylabel(r"$p_\mathrm{hist}$ (2a)")
ax.invert_yaxis()                            # AA plateau increases upward
ax.set_title(r"c  $\Delta D$, age × plateau" "\n" "(rise timing re-optimised)")
ax = fig.add_subplot(gs[1, 2:])
xs, ys, im = heat(ax, 3, 2, "")
ax.set_xticks(range(len(xs)), [f"{x / 1e3:.0f}" for x in xs]); ax.set_yticks(range(len(ys)), [f"{y / 1e3:.0f}" for y in ys])
ax.set_xlabel(r"$t_\mathrm{rise}$ (ky)"); ax.set_ylabel(r"$t_\mathrm{end}$ (ka)")
ax.set_title(r"d  $\Delta D$, rise duration × end of rise" "\n" "(age and plateau re-optimised)")
cb = fig.colorbar(im, ax=fig.axes[-2:], fraction=0.025, pad=0.02)
cb.set_label(r"$\Delta D$ from the best fit (darker is better; – not simulated)")

# ---------------- (e) profiles with calibrated cutoffs ----------------
lab = {"t_inv": (r"$t_\mathrm{inv}$ (ky)", 1e3), "p_hist": (r"$p_\mathrm{hist}$ (2a)", 1),
       "t_decline": (r"$t_\mathrm{end}$ (ka)", 1e3), "t_fall": (r"$t_\mathrm{rise}$ (ky)", 1e3)}
for j, p in enumerate(PARS):
    ax = fig.add_subplot(gs[2, j])
    v = np.array(sorted(prof[p])); d = np.array([prof[p][x] for x in v])
    if p == "p_hist":                       # plot the AA plateau, 1 - p_BB
        v = 1 - v
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
    if j == 1:
        ax.legend(loc="upper center", fontsize=7)
    if j == 0:
        ax.set_ylabel(r"profile $\Delta D$ (+0.1, log)")
        ax.set_title(r"e  Profiles, shaded below the bootstrap-calibrated cutoff (dashed), dotted line at $\Delta D=2$",
                     loc="left")
        pass

# ---------------- (f) why the neighbours are rejected ----------------
ax = fig.add_subplot(gs[3, :])
rows = [(f"best fit ({fmt(cell(*best))})", cell(*best)),
        ("barrier 775 ky", cell(775e3, None, None, None)),
        ("barrier 875 ky", cell(875e3, None, None, None)),
        ("barrier 700 ky", cell(700e3, None, None, None)),
        ("plateau 0.40", cell(None, 0.60, None, None)),
        ("plateau 0.18", cell(None, 0.82, None, None)),
        ("rise took 400 ky", cell(None, None, None, 400e3)),
        ("rise ended 50 ka", cell(None, None, 50e3, None)),
        ("rise ends today", cell(None, None, 0.0, None))]
y = np.arange(len(rows))[::-1]
for yi, (name, k) in zip(y, rows):
    left = 0.0
    for s_, col in enumerate(STAT_COL):
        w = z2[k, s_]
        ax.barh(yi, w, left=left, height=0.62, color=col, edgecolor="white", lw=1.5)
        if w > 1.2:
            ax.text(left + w / 2, yi, f"{w:.1f}", ha="center", va="center", fontsize=7.5, color=C["ink"])
        left += w
    ax.text(left + 0.6, yi, r"$D$ " + f"{left:.1f}" + ("" if k == cell(*best) else f"   ({fmt(k)})"), va="center",
            fontsize=8, color=C["ink"])
ax.set_yticks(y, [r[0] for r in rows])
ax.set_xlabel(r"sequential contribution to $D$ (each statistic given those before it), other parameters re-optimised")
ax.set_xlim(0, max(chi[r[1]] for r in rows) * 1.45)
ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(STAT_COL, STAT_LAB_F)], loc="lower right",
          ncol=3, fontsize=8)
ax.set_title("f  Which statistic rules out each neighbouring value", loc="left")
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)

fig.savefig(OUT, dpi=170, bbox_inches="tight")
print("wrote", OUT)
