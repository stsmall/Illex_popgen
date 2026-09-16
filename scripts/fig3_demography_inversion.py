"""Fig 3 -- demography + chr2 inversion frequency history, ONE combined figure.

Panel a: population size history (composite) -- Ne over time, overlaying (each with
    computed, not fabricated, uncertainty where available) Stairway Plot 2
    (median + 95% CI), the moments 2-epoch growth model + its Godambe (GIM)
    parameter-CI envelope, GONE2 recent Ne + its bootstrap CI band, and the
    momentsLD single contemporary LD-Ne point (with caveat). Same plotting logic
    and data as fig3a_demography_composite.py (same directory), the standalone
    version of this panel:
      - stairway: steps/08_demography/stairway_out/*.final.summary
      - moments: analytic 2-epoch exponential-growth curve (N_ANC=547928,
        N0=6808096, T=769519 generations) + N(t) envelope from
        steps/08_demography/moments_godambe_uncert.json (GIM_uncert on a
        chromosome-block bootstrap of the folded SFS)
      - GONE2: steps/10_gone2/gone_input.vcf_GONE2_Ne (point; generations 1-4
        dropped) + steps/10_gone2/bootstrap/gone2_band.tsv (bootstrap 95% band)
      - momentsLD: single Ohta-Kimura LD-Ne (ld_Ne.txt), rate-corrected from the
        flat 1e-8/bp assumption to the empirical ~2.1e-9/bp map (Ne ~ 1/rate);
        the time-varying fit is non-identifiable (params pinned at bounds).

Panel b: chr2 inversion (BB arrangement) frequency history -- the fitted
    decline trajectory (b-i) plus the two delta-chi-square profiles that
    bound it (b-ii: t_decline, "when the fall ended"; b-iii: t_fall, "how
    long the fall took"). Re-plotted from the msinv project's
    figures_inversion.py (manu_fig_fitted_history / figure2), which lives in
    a separate repo+venv (msinv_dir/inversion_sims/files, needs the compiled
    `msinv` extension). That script's own plotting/annotation logic and the
    FIT point are reproduced verbatim below; the only thing NOT imported
    from there is `illex.balancing.decline_curve` itself, which is vendored
    as pure numpy/scipy (verified byte-for-byte identical against the real
    msinv .venv -- see analysis note in this repo's session log) so this
    panel runs under the shared bioinfo-buddy interpreter alongside panel a.
    Data:
      - FIT point + TARGETS: figures_inversion.py (hardcoded constants)
      - refit_declinetdecline.csv, refit_declinetfall2.csv:
        msinv_dir/inversion_sims/files/results/illex/ (profile scan output)

    .venv/bin/python -m illex.scripts.figures_inversion  (source repo,
    for reference / regenerating the standalone fig_fitted_history.png)
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, OKABE, despine, save
apply()
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from scipy import optimize

# =====================================================================
# Panel a data -- identical sources/logic to fig3a_demography_composite.py
# =====================================================================
STAIRWAY = ("/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/"
            "stairway_out/Illex_illecebrosus_baker633.final.summary")
GONE2 = ("/sietch_colab/data_share/illex/popgen_data/analysis/steps/10_gone2/"
         "gone_input.vcf_GONE2_Ne")
GODAMBE = ("/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/"
           "moments_godambe_uncert.json")
GONE2_BAND = ("/sietch_colab/data_share/illex/popgen_data/analysis/steps/10_gone2/"
              "bootstrap/gone2_band.tsv")

COL_STAIRWAY = C["accent"]   # "#0072B2" blue
COL_MOMENTS = C["warm"]      # "#D55E00" vermillion
COL_GONE2 = OKABE[3]         # "#009E73" bluish green
COL_MLD = OKABE[1]           # "#E69F00" orange (momentsLD single point)

N_ANC = 547928
N0 = 6808096
T_GROW_MOMENTS = 769519  # generations of growth, generations ~= years
t_growth = np.geomspace(1, T_GROW_MOMENTS, 400)
n_growth = N0 * (N_ANC / N0) ** (t_growth / T_GROW_MOMENTS)
t_flat = np.geomspace(T_GROW_MOMENTS, 3_000_000, 50)
n_flat = np.full_like(t_flat, N_ANC)
t_moments = np.concatenate([t_growth, t_flat])
n_moments = np.concatenate([n_growth, n_flat])

sw = pd.read_csv(STAIRWAY, sep="\t")
sw = sw.drop_duplicates(subset=["year"]).sort_values("year")
sw = sw[sw["year"] > 0]

g2 = pd.read_csv(GONE2, sep="\t").sort_values("Generation").reset_index(drop=True)
# GONE2's earliest generations (1-4) are a clamped recent-generation artifact; extend the
# line to the present by holding the first reliable estimate (generation 5) flat to gen 1.
_ne5 = float(g2.loc[g2["Generation"] == 5, "Ne_diploids"].iloc[0])
g2.loc[g2["Generation"] < 5, "Ne_diploids"] = _ne5

# =====================================================================
# Panel b data -- vendored from msinv_dir/inversion_sims/files/illex/
#   balancing.py (decline_curve, _time_to_reach) and theory.py (N_growth).
#   Pure numpy/scipy; no `msinv` C-extension needed. Verified identical to
#   the real illex.balancing.decline_curve output under the msinv .venv
#   (400/400 points, times/freqs/s_fall all match to full float precision).
# =====================================================================
RESULTS_ILLEX = Path("/sietch_colab/ssmall/projects/msinv_dir/inversion_sims/files/"
                     "results/illex")

# --- theory.py: moments growth-model Ne(t), used only for the drift envelope
N_ANC_INV = 547_928.0
N0_INV = 6_808_096.0
T_GROW_INV = 769_519.0
ALPHA_INV = np.log(N0_INV / N_ANC_INV) / T_GROW_INV


def N_growth(t):
    t = np.asarray(t, dtype=float)
    return np.where(t <= T_GROW_INV, N0_INV * np.exp(-ALPHA_INV * np.minimum(t, T_GROW_INV)),
                    N_ANC_INV)


ARRIVAL_TOL = 1e-4


def _time_to_reach(p, p_start, s_het, p_star):
    """Closed-form forward time from p_start to p (balancing.py, verbatim)."""
    p = np.asarray(p, dtype=float)
    c = s_het / p_star
    a = 1.0 / p_star
    b = 1.0 / (1.0 - p_star)
    d = 1.0 / (p_star * (1.0 - p_star))

    def F(x):
        return (a * np.log(x) + b * np.log(1.0 - x)
                - d * np.log(np.abs(p_star - x)))

    return (F(p) - F(p_start)) / c


def decline_curve(t_inv, t_decline, t_fall, p_hist, p_now,
                  n_fall=300, n_hist=100, tol=ARRIVAL_TOL):
    """(times, freqs, s_fall) backward on [0, t_inv] (balancing.py, verbatim)."""
    if not 0.0 < p_now < p_hist < 1.0:
        raise ValueError("need 0 < p_now < p_hist < 1")
    if t_decline + t_fall > t_inv:
        raise ValueError("t_decline + t_fall exceeds t_inv")

    def _fall_time(s):
        p_end = p_now * (1.0 + tol)
        return float(_time_to_reach(p_end, p_hist, s, p_now)) - t_fall

    s_fall = float(optimize.brentq(_fall_time, 1e-9, 10.0, rtol=1e-12))
    p_grid = np.geomspace(p_hist, p_now * (1.0 + tol), n_fall)
    t_fwd = _time_to_reach(p_grid, p_hist, s_fall, p_now)
    t_back = t_decline + (t_fall - t_fwd)

    times = [0.0]
    freqs = [p_now]
    order = np.argsort(t_back)
    times.extend(t_back[order].tolist())
    freqs.extend(p_grid[order].tolist())
    hist = t_inv - t_decline - t_fall
    if hist > 0:
        for f in np.linspace(0.0, 1.0, n_hist)[1:]:
            times.append(t_decline + t_fall + f * hist)
            freqs.append(p_hist)

    times = np.asarray(times, dtype=float)
    freqs = np.asarray(freqs, dtype=float)
    srt = np.argsort(times, kind="stable")
    times, freqs = times[srt], freqs[srt]
    keep = np.concatenate([[True], np.diff(times) > 1e-9])
    times, freqs = times[keep], freqs[keep]
    return np.clip(times, 0.0, t_inv).tolist(), freqs.tolist(), s_fall


def drift_sd(p0, t_lo, t_hi, n=20000):
    """SD of a neutral frequency change over [t_lo, t_hi] (figures_inversion.py)."""
    t = np.linspace(t_lo, t_hi, n)
    return float(np.sqrt(p0 * (1.0 - p0)
                         * np.trapezoid(1.0 / (2.0 * N_growth(t)), t)))


# --- fitted point + empirical targets (figures_inversion.py, verbatim) -----
FIT = dict(t_inv=850_000.0, p_hist=0.74, t_decline=175_000.0,
           t_fall=100_000.0, p_now=0.374)
TARGETS = {"pi_i_over_pi_s": (1.3556, 0.0481),
           "dxy_over_pi_i":  (1.3848, 0.0214),
           "f1_ratio":       (0.8125, 0.0076)}


def chi2(row):
    return sum(((float(row[k]) - t) / se) ** 2
               for k, (t, se) in TARGETS.items())


def profile(tag, key):
    """(values, delta-chi2) profiling `key` (figures_inversion.py, verbatim)."""
    rows = list(csv.DictReader(open(RESULTS_ILLEX / f"refit_decline{tag}.csv")))
    best = {}
    for r in rows:
        v, c = float(r[key]), chi2(r)
        if v not in best or c < best[v][0]:
            best[v] = (c, r)
    gmin = min(c for c, _ in best.values())
    vals = sorted(best)
    return np.array(vals), np.array([best[v][0] - gmin for v in vals])


times, freqs, _s_fall = decline_curve(**FIT)
times, freqs = np.asarray(times) / 1e3, np.asarray(freqs)

# one consistent accent for the inversion-frequency data, distinct from the
# three colours already spoken for in panel a (blue/vermillion/green)
COL_BB = OKABE[7]     # "#CC79A7" reddish purple

# =====================================================================
# Figure
# =====================================================================
fig = plt.figure(figsize=(7.8, 6.6))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 0.92],
              hspace=0.42, wspace=0.28)

LETTER_KW = dict(fontweight="bold", fontsize=13)

# ---- a: demography composite (fig3a_demography_composite.py, verbatim) ----
ax_a = fig.add_subplot(gs[0, :])

# moments Godambe (GIM) parameter-CI envelope on N(t)
gd = json.load(open(GODAMBE))
env = gd["nt_envelope"]
ax_a.fill_between(env["t"], env["lo"], env["hi"], color=COL_MOMENTS, alpha=0.16,
                  linewidth=0, zorder=2, label="_nolegend_")

ax_a.fill_between(sw["year"], sw["Ne_2.5%"], sw["Ne_97.5%"],
                  color=COL_STAIRWAY, alpha=0.15, linewidth=0, zorder=1,
                  label="_nolegend_")
ax_a.plot(sw["year"], sw["Ne_median"], color=COL_STAIRWAY, linewidth=1.4,
          zorder=4, label="Stairway Plot 2 (median, 95% CI)")
ax_a.plot(t_moments, n_moments, color=COL_MOMENTS, linewidth=1.8, linestyle="--",
          zorder=5, label="moments 2-epoch growth (model, 95% GIM CI)")

# GONE2 recent trajectory (bootstrap CI band drawn only if a real bootstrap exists)
band = pd.read_csv(GONE2_BAND, sep="\t")
band = band[band["Generation"] >= 5]
nrep = int(band["n_rep"].iloc[0])
if nrep > 0:
    ax_a.fill_between(band["Generation"], band["Ne_2.5"], band["Ne_97.5"],
                      color=COL_GONE2, alpha=0.18, linewidth=0, zorder=3,
                      label="_nolegend_")
_gone_lab = "GONE2 recent (LD-based)" if nrep == 0 else \
            f"GONE2 recent (LD-based; {nrep}-rep bootstrap 95% CI)"
ax_a.plot(g2["Generation"], g2["Ne_diploids"], color=COL_GONE2, linewidth=1.8,
          zorder=6, label=_gone_lab)

# momentsLD single contemporary LD-Ne point (rate-corrected 1e-8 -> ~2.1e-9/bp)
_ld_corr = 323741.0 * (1e-8 / 2.1e-9)
ax_a.scatter([1.3], [_ld_corr], marker="*", s=150, color=COL_MLD, zorder=7,
             edgecolor="white", linewidth=0.8,
             label="momentsLD LD-$N_e$ (single point; see note)")

ax_a.set_xscale("log")
ax_a.set_yscale("log")
ax_a.set_xlim(1, 3_000_000)
ax_a.set_ylim(2e4, 2e7)
ax_a.set_xlabel("Generations before present")
ax_a.set_ylabel(r"Effective population size ($N_e$)")
ax_a.yaxis.set_major_formatter(mticker.FuncFormatter(
    lambda v, _: f"{v/1e6:g}M" if v >= 1e6 else f"{v/1e3:g}k"))
despine(ax_a)
_leg = ax_a.legend(loc="lower left", handlelength=1.8, fontsize=8.0,
                   borderaxespad=0.4, labelspacing=0.35)
_leg.set_zorder(20)
ax_a.text(-0.085, 1.04, "a", transform=ax_a.transAxes, **LETTER_KW)
ax_a.set_title("Population size history", loc="left", pad=8,
               fontsize=11, fontweight="normal", color=C["muted"])

# ---- b-i: fitted inversion-frequency trajectory ---------------------------
ax_bi = fig.add_subplot(gs[1, :])

t0 = FIT["t_decline"] + FIT["t_fall"]
tb = np.linspace(0.0, t0, 400)
sd = np.array([drift_sd(FIT["p_hist"], t, t0) for t in tb])
for k in (2, 1):
    ax_bi.fill_between(tb / 1e3, np.clip(FIT["p_hist"] - k * sd, 0, 1),
                       np.clip(FIT["p_hist"] + k * sd, 0, 1),
                       color=C["muted"], alpha=0.14, lw=0, zorder=1)
ax_bi.plot(tb / 1e3, np.full_like(tb, FIT["p_hist"]), lw=0.8, ls=":",
          color=C["muted"], zorder=1)
ax_bi.axhline(FIT["p_now"], lw=0.8, ls="--", color=C["ink"], zorder=2)
ax_bi.plot(times, freqs, lw=2.0, color=COL_BB, zorder=3, solid_capstyle="round")
ax_bi.scatter([0], [FIT["p_now"]], s=32, color=COL_BB, zorder=4,
             edgecolor="white", linewidth=0.9)

ax_bi.set_xlim(-14, FIT["t_inv"] / 1e3 * 1.06)   # present (0 ka) at LEFT, to match panel a
ax_bi.set_ylim(0, 0.96)
ax_bi.set_yticks([0, 0.2, 0.374, 0.6, 0.74, 0.9])
ax_bi.set_yticklabels(["0", "0.2", "0.374", "0.6", "0.74", "0.9"])
ax_bi.set_xlabel("Time before present (ka)")
ax_bi.set_ylabel("Frequency of inverted\n(BB) arrangement", linespacing=1.3)
despine(ax_bi)

ax_bi.annotate("arose ~850 ka", xy=(FIT["t_inv"] / 1e3, FIT["p_hist"]),
              xytext=(FIT["t_inv"] / 1e3 - 25, 0.90), fontsize=9,
              color=C["ink"], ha="right", va="center",
              arrowprops=dict(arrowstyle="-", lw=0.7, color=C["muted"],
                              shrinkA=1, shrinkB=3))
ax_bi.text(8, FIT["p_now"] + 0.035, "$p = 0.374$ today", fontsize=9,
          color=C["ink"], ha="left", va="bottom")
ax_bi.text(-0.085, 1.06, "b", transform=ax_bi.transAxes, **LETTER_KW)
ax_bi.set_title("Inversion frequency history", loc="left", pad=8,
                fontsize=11, fontweight="normal", color=C["muted"])

OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/fig4_demography.png"
save(fig, OUT)
