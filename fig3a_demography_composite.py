"""Fig 3A -- composite demographic-history figure (Ne over time), with uncertainty.

Overlays, each with principled (computed, not fabricated) uncertainty where available:
  - Stairway Plot 2: median Ne(t) + 2.5-97.5% CI band, from
    steps/08_demography/stairway_out/Illex_illecebrosus_baker633.final.summary
  - moments 2-epoch exponential-growth model: the analytic model curve PLUS a shaded
    parameter-CI envelope. The envelope is the 2.5-97.5% band of N(t) implied by the
    Godambe Information Matrix (GIM) uncertainty on (nu, T, theta), computed by
    moments.Godambe.GIM_uncert on a chromosome-block bootstrap of the folded SFS
    (steps/08_demography/moments_godambe.py -> moments_godambe_uncert.json).
  - GONE2 recent-Ne: median trajectory PLUS a 2.5-97.5% bootstrap CI band, from
    re-running GONE2 on chromosome-block-bootstrap replicates
    (steps/10_gone2/bootstrap/gone2_band.tsv; point line = the original single run,
    steps/10_gone2/gone_input.vcf_GONE2_Ne). Generations 1-4 dropped (known GONE2
    recent-edge artifact).
  - momentsLD: a single contemporary LD-effective-Ne point at the present, shown with
    an explicit caveat. The momentsLD *time-varying* fit is non-identifiable here
    (optimizer pins parameters against their bounds; steps/09_momentsld/momentsld_fit.json),
    so only the single Ohta-Kimura LD-Ne is shown -- rate-corrected from the flat
    1e-8/bp assumption (ld_Ne.txt) to the empirical ~2.1e-9/bp map (Ne ~ 1/rate).

x = time before present in generations (~= years), log scale.
y = effective population size Ne, log scale.
"""
import json
import sys

sys.path.insert(0, "/sietch_colab/data_share/illex/popgen_data/analysis/manuscript")
from figstyle import apply, C, OKABE, despine, save
apply()
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

D8 = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography"
D10 = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/10_gone2"
STAIRWAY = f"{D8}/stairway_out/Illex_illecebrosus_baker633.final.summary"
GODAMBE = f"{D8}/moments_godambe_uncert.json"
GONE2_POINT = f"{D10}/gone_input.vcf_GONE2_Ne"
GONE2_BAND = f"{D10}/bootstrap/gone2_band.tsv"
OUT = ("/sietch_colab/data_share/illex/popgen_data/analysis/manuscript/figures/"
       "fig_demography_composite.png")

# fixed Okabe-Ito colour assignment, one per method (not a rainbow)
COL_STAIRWAY = C["accent"]   # "#0072B2" blue
COL_MOMENTS = C["warm"]      # "#D55E00" vermillion
COL_GONE2 = OKABE[3]         # "#009E73" bluish green
COL_MLD = OKABE[1]           # "#E69F00" orange (momentsLD point)


def draw_panel_a(ax, legend_loc="lower left"):
    # ---- moments 2-epoch exponential growth model (analytic point curve) ----
    N_ANC, N0, T = 547928, 6808096, 769519
    t_growth = np.geomspace(1, T, 400)
    n_growth = N0 * (N_ANC / N0) ** (t_growth / T)
    t_flat = np.geomspace(T, 3_000_000, 50)
    t_moments = np.concatenate([t_growth, t_flat])
    n_moments = np.concatenate([n_growth, np.full_like(t_flat, N_ANC)])

    # ---- moments Godambe parameter-CI envelope ----
    gd = json.load(open(GODAMBE))
    env = gd["nt_envelope"]
    ax.fill_between(env["t"], env["lo"], env["hi"], color=COL_MOMENTS, alpha=0.16,
                    linewidth=0, zorder=2, label="_nolegend_")

    # ---- Stairway Plot 2 ----
    sw = pd.read_csv(STAIRWAY, sep="\t").drop_duplicates(subset=["year"]).sort_values("year")
    sw = sw[sw["year"] > 0]
    ax.fill_between(sw["year"], sw["Ne_2.5%"], sw["Ne_97.5%"], color=COL_STAIRWAY,
                    alpha=0.15, linewidth=0, zorder=1, label="_nolegend_")
    ax.plot(sw["year"], sw["Ne_median"], color=COL_STAIRWAY, linewidth=1.4, zorder=4,
            label="Stairway Plot 2 (median, 95% CI)")

    # moments point curve on top of its envelope
    ax.plot(t_moments, n_moments, color=COL_MOMENTS, linewidth=1.8, linestyle="--",
            zorder=5, label="moments 2-epoch growth (model, 95% GIM CI)")

    # ---- GONE2 median + bootstrap CI band ----
    band = pd.read_csv(GONE2_BAND, sep="\t")
    band = band[band["Generation"] >= 5]
    ax.fill_between(band["Generation"], band["Ne_2.5"], band["Ne_97.5"], color=COL_GONE2,
                    alpha=0.18, linewidth=0, zorder=3, label="_nolegend_")
    g2 = pd.read_csv(GONE2_POINT, sep="\t")
    g2 = g2[g2["Generation"] >= 5]
    nrep = int(band["n_rep"].iloc[0])
    ax.plot(g2["Generation"], g2["Ne_diploids"], color=COL_GONE2, linewidth=1.8, zorder=6,
            label=f"GONE2 recent (LD-based; {nrep}-rep bootstrap 95% CI)")

    # ---- momentsLD single contemporary LD-Ne point (rate-corrected) ----
    ld_flat = 323741.0
    flat_rate, real_rate = 1e-8, 2.1e-9
    ld_corr = ld_flat * (flat_rate / real_rate)   # Ne ~ 1/rate  -> ~1.54e6
    ax.scatter([1.3], [ld_corr], marker="*", s=150, color=COL_MLD, zorder=7,
               edgecolor="white", linewidth=0.8,
               label="momentsLD LD-$N_e$ (single point; see note)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, 3_000_000)
    ax.set_ylim(2e4, 2e7)
    ax.set_xlabel("Generations before present")
    ax.set_ylabel(r"Effective population size ($N_e$)")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f"{v/1e6:g}M" if v >= 1e6 else f"{v/1e3:g}k"))
    despine(ax)
    leg = ax.legend(loc=legend_loc, handlelength=1.8, fontsize=8.2,
                    borderaxespad=0.4, labelspacing=0.35)
    leg.set_zorder(20)
    # caveat note for the momentsLD point (non-identifiable time-varying fit)
    ax.annotate("momentsLD time-varying fit non-identifiable\n"
                "(params pinned at bounds); only single LD-$N_e$ shown,\n"
                "rate-corrected to ~2.1e-9/bp map",
                xy=(1.3, ld_corr), xytext=(28, 2.7e5), fontsize=6.6, color=C["muted"],
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C["muted"],
                                shrinkA=2, shrinkB=4))
    return ld_corr


if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    draw_panel_a(ax)
    save(fig, OUT)
