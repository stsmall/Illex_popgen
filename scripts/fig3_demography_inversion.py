"""Fig 3 -- demography + chr2 inversion frequency history, ONE combined figure.

Panel a: population size history (composite) -- Ne over time, overlaying (each with
    computed, not fabricated, uncertainty where available) Stairway Plot 2
    (median + 95% CI), the moments 2-epoch growth model + its Godambe (GIM)
    parameter-CI envelope, and GONE2 recent Ne + its bootstrap CI band. Same plotting logic
    and data as fig3a_demography_composite.py (same directory), the standalone
    version of this panel:
      - stairway: steps/08_demography/stairway_out/*.final.summary
      - moments: analytic 2-epoch exponential-growth curve (N_ANC=547928,
        N0=6808096, T=769519 generations) + N(t) envelope from
        steps/08_demography/moments_godambe_uncert.json (GIM_uncert on a
        chromosome-block bootstrap of the folded SFS)
      - GONE2: steps/10_gone2/gone_input.vcf_GONE2_Ne (all 43 autosomes, -r 1 i.e. an assumed
        uniform 1 cM/Mb; generations 1-4 held at the generation-5 value) + its bootstrap band if one
        exists. A rerun at the ReLERNN rate (-r 0.24, steps/10_gone2/r024/) is not shown: GONE2 needs
        every chromosome to span >20 cM, so at 0.24 cM/Mb only 11 autosomes qualify and that run is
        unstable in recent generations.

Panel b: chr2 inversion frequency history of the derived 2a arrangement under the balancing model (B) with
    gene flux, the model favoured by the ABC-SMC model choice, with parameters from neural posterior estimation
    (NPE). Pointwise quantiles (50% and 95%) and median of 400 trajectories re-simulated from NPE posterior draws
    on the inversion's neutral baseline. Individual trajectories are not drawn (they are in Supplementary Fig. S14,
    scripts/figS14c_abc_trajectory.py). Data: data/inversion_abc/traj_quantiles_npe_f1.csv, written by
    illex/abc/npe/plot_npe.py in the msinv repository (branch illex/trajectory-abc).
"""
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
# Panel b data -- NPE posterior trajectories (model B, gene flux)
# =====================================================================
TRAJ = Path(__file__).resolve().parents[1] / "data" / "inversion_abc" / "traj_quantiles_npe_f1.csv"
tq = pd.read_csv(TRAJ)
P_NOW = 0.618  # 2a frequency over 621 karyotyped individuals (242/284/95)

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
nrep = 0
if GONE2_BAND is not None:
    band = pd.read_csv(GONE2_BAND, sep="\t")
    band = band[band["Generation"] >= 5]
    nrep = int(band["n_rep"].iloc[0])
if nrep > 0:
    ax_a.fill_between(band["Generation"], band["Ne_2.5"], band["Ne_97.5"],
                      color=COL_GONE2, alpha=0.18, linewidth=0, zorder=3,
                      label="_nolegend_")
_gone_lab = "GONE2 recent (LD-based, 1 cM/Mb)" if nrep == 0 else \
            f"GONE2 recent (LD-based; {nrep}-rep bootstrap 95% CI)"
ax_a.plot(g2["Generation"], g2["Ne_diploids"], color=COL_GONE2, linewidth=1.8,
          zorder=6, label=_gone_lab)



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

# ---- b: posterior trajectory of the derived arrangement -------------------
ax_bi = fig.add_subplot(gs[1, :])
tk = tq["t_gen"] / 1e3
ax_bi.fill_between(tk, tq["q025"], tq["q975"], color=COL_BB, alpha=0.16, lw=0, zorder=1, label="95%")
ax_bi.fill_between(tk, tq["q25"], tq["q75"], color=COL_BB, alpha=0.32, lw=0, zorder=2, label="50%")
ax_bi.plot(tk, tq["q50"], lw=2.0, color=COL_BB, zorder=4, solid_capstyle="round", label="median")
ax_bi.axhline(P_NOW, lw=0.8, ls="--", color=C["ink"], zorder=2)
ax_bi.set_xlim(-14, 1000)   # present (0 ka) at LEFT, to match panel a
ax_bi.set_ylim(0, 0.8)
ax_bi.set_yticks([0, 0.2, 0.4, P_NOW])
ax_bi.set_yticklabels(["0", "0.2", "0.4", f"{P_NOW:.3f}"])
ax_bi.set_xlabel("Time before present (ka)")
ax_bi.set_ylabel("Frequency of derived\n(2a) arrangement", linespacing=1.3)
despine(ax_bi)
ax_bi.legend(loc="center right", fontsize=8.0, title="posterior", title_fontsize=8.0)
ax_bi.text(8, P_NOW + 0.035, f"$p = {P_NOW:.3f}$ today", fontsize=9, color=C["ink"], ha="left", va="bottom")
ax_bi.text(-0.085, 1.06, "b", transform=ax_bi.transAxes, **LETTER_KW)
ax_bi.set_title("Inversion frequency history", loc="left", pad=8,
                fontsize=11, fontweight="normal", color=C["muted"])

OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/fig4_demography.png"
save(fig, OUT)
