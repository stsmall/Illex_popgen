#!/usr/bin/env python
"""Godambe (GIM) parameter uncertainty for the adopted moments `growth` model.

Block bootstrap = resample CHROMOSOMES with replacement (43 autosomal blocks, excl.
chr2 inversion + chr42 sex), re-tally the folded genome-wide SFS, and feed the
bootstrap spectra to moments.Godambe.GIM_uncert. Chromosomes are the natural
LD-respecting block: sites within a chromosome are linked, so the appropriate
bootstrap chunk is chromosome-scale, exactly as the moments/dadi Godambe docs
recommend. Re-running the ANGSD/realSFS GL EM per replicate is infeasible here
(1266 haploids, ~1e9 sites), so the already-estimated per-chromosome folded SFS
point estimates (thetas/baker.<c>.sfs) are resampled with replacement and re-summed.

To keep moments' SFS integration fast and numerically stable, both the observed and
the bootstrap folded spectra are PROJECTED down to n=200 haploids. The demographic
parameters (nu, T) and the population-scaled theta are sample-size independent, so
projection changes neither the inference target nor Nref; it only makes each
`growth()` evaluation ~100x cheaper (0.03 s vs seconds at n=1266) and avoids the
optimizer/finite-difference steps stalling in stiff-integration regions.

Point estimate = the adopted full-data `growth` fit (moments_fit_results.tsv:
N_ANC=547928 -> N0=6808096 over T=769519 gen). The GIM supplies the (relative)
parameter covariance, which is applied about that point and Monte-Carlo propagated to
CIs on N_ANC, N0 (present Ne), T (generations) and the expansion factor.
 -> moments_godambe_uncert.json
"""
import json
import numpy as np
import moments

D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography"
TH = f"{D}/thetas"
MU = 3e-9
GEN = 1.0
L = 1340915910                      # accessible sites with data (autosomes excl chr2,chr42)
N = 1266                            # haploid sample size (633 diploids)
PROJ = 200                          # projected haploid size for fast/stable moments
POOL = [c for c in range(1, 46) if c not in (2, 42)]   # 43 autosomal blocks
NBOOT = 100
SEED = 20240914
rng = np.random.default_rng(SEED)

# adopted full-data growth fit (moments_fit_results.tsv, row all_bins/growth)
N_ANC_PUB = 547928.0
N0_PUB = 6808096.0
T_YEARS_PUB = 769519.0
NU_PUB = N0_PUB / N_ANC_PUB                      # 12.4251
T_PUB = T_YEARS_PUB / (2 * N_ANC_PUB)           # 0.70223 (coalescent units)
THETA_PUB = 4 * MU * L * N_ANC_PUB

# ---- load per-chromosome folded SFS (realSFS -fold output: 1267 floats each) ----
chrom_sfs = {}
for c in POOL:
    with open(f"{TH}/baker.{c}.sfs") as fh:
        v = np.array(fh.read().split(), dtype=float)
    assert v.size == 1267, f"chr{c}: {v.size} entries"
    chrom_sfs[c] = v
chroms = list(chrom_sfs)
print(f"loaded {len(chroms)} per-chromosome SFS blocks")


def make_fs(folded_counts):
    """Full folded moments.Spectrum (n=1266), then project to PROJ."""
    data = np.zeros(N + 1)
    data[:len(folded_counts)] = folded_counts
    fs = moments.Spectrum(data, mask_corners=False)
    fs.folded = True
    fs.mask[0] = True
    fs.mask[N // 2 + 1:] = True
    fs.mask[N] = True
    return fs.project([PROJ])


# ---- observed genome-wide folded SFS (bins 0..633) -> projected ----
tot = np.sum([chrom_sfs[c] for c in chroms], axis=0)
fs = make_fs(tot[:N // 2 + 1].copy())
ns = fs.sample_sizes


def growth_folded(p, ns):
    return moments.Demographics1D.growth(p, ns).fold()


# ---- polish the fit on the projected observed SFS (fast, stable) ----
lb, ub = [1e-3, 1e-4], [1e4, 10]
popt = moments.Inference.optimize_log([NU_PUB, T_PUB], fs, growth_folded,
                                      lower_bound=lb, upper_bound=ub, maxiter=100, verbose=0)
model = growth_folded(popt, ns)
ll_opt = moments.Inference.ll_multinom(model, fs)
theta_proj = moments.Inference.optimal_sfs_scaling(model, fs)
Nref_proj = theta_proj / (4 * MU * L)
print(f"projected fit: nu={popt[0]:.4f} T={popt[1]:.5f} theta={theta_proj:.4g} "
      f"logL={ll_opt:.1f} Nref={Nref_proj:,.0f}")
print(f"adopted point: N_ANC={N_ANC_PUB:,.0f} N0={N0_PUB:,.0f} T_gen={T_YEARS_PUB:,.0f} "
      f"expansion={NU_PUB:.2f}")

# ---- bootstrap SFS replicates: resample chromosomes with replacement -> project ----
boots = []
for b in range(NBOOT):
    pick = rng.choice(chroms, size=len(chroms), replace=True)
    bsum = np.sum([chrom_sfs[c] for c in pick], axis=0)
    boots.append(make_fs(bsum[:N // 2 + 1].copy()))
print(f"built {len(boots)} projected bootstrap SFS replicates (chromosome block bootstrap)")

# ---- Godambe Information Matrix uncertainty (at the projected optimum) ----
uncerts_lin = moments.Godambe.GIM_uncert(
    growth_folded, boots, popt, fs, log=False, multinom=True)
uncerts_log, GIM_log = moments.Godambe.GIM_uncert(
    growth_folded, boots, popt, fs, log=True, multinom=True, return_GIM=True)
cov_log = np.linalg.inv(GIM_log)               # covariance of [log nu, log T, log theta]
print("GIM SDs (linear):  nu=%.4g  T=%.4g  theta=%.4g" % tuple(uncerts_lin))
print("GIM SDs (log/rel): nu=%.4g  T=%.4g  theta=%.4g" % tuple(uncerts_log))

# ---- Monte-Carlo propagation, centred on the adopted point, with GIM log-cov ----
mean_log = np.log([NU_PUB, T_PUB, THETA_PUB])
NDRAW = 200000
draws_log = rng.multivariate_normal(mean_log, cov_log, size=NDRAW)
draws = np.exp(draws_log)
nu_d, T_d, th_d = draws[:, 0], draws[:, 1], draws[:, 2]
Nref_d = th_d / (4 * MU * L)
derived = {
    "nu": nu_d, "T_coal": T_d, "theta": th_d,
    "N_ANC": Nref_d, "N0_present": nu_d * Nref_d,
    "T_generations": 2 * Nref_d * T_d * GEN, "expansion_factor": nu_d,
}
point_vals = {
    "nu": NU_PUB, "T_coal": T_PUB, "theta": THETA_PUB,
    "N_ANC": N_ANC_PUB, "N0_present": N0_PUB,
    "T_generations": T_YEARS_PUB, "expansion_factor": NU_PUB,
}


def summ(x, point):
    return dict(point=float(point), mean=float(np.mean(x)), sd=float(np.std(x, ddof=1)),
                ci95=[float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))],
                cv=float(np.std(x, ddof=1) / point) if point else None)


# ---- N(t) parameter-CI envelope for the figure ----
tgrid = np.geomspace(1.0, 3_000_000.0, 300)
sub = draws[rng.choice(draws.shape[0], size=5000, replace=False)]
curves = np.empty((sub.shape[0], tgrid.size))
for i, (nu_i, T_i, th_i) in enumerate(sub):
    Nref_i = th_i / (4 * MU * L)
    N0_i = nu_i * Nref_i
    Tg_i = 2 * Nref_i * T_i * GEN
    curves[i] = np.where(tgrid <= Tg_i, N0_i * (Nref_i / N0_i) ** (tgrid / Tg_i), Nref_i)
env_lo = np.percentile(curves, 2.5, axis=0)
env_hi = np.percentile(curves, 97.5, axis=0)
env_med = np.percentile(curves, 50, axis=0)

results = {
    "method": "moments.Godambe.GIM_uncert on chromosome-block-bootstrap folded SFS (projected to n=%d)" % PROJ,
    "model": "growth (2-epoch exponential)",
    "n_bootstrap": NBOOT,
    "n_blocks": len(chroms),
    "block_type": "chromosome (43 autosomes, excl chr2 inversion + chr42 sex)",
    "projected_n_hap": PROJ, "full_n_hap": N,
    "mu": MU, "gen_time": GEN, "L": L,
    "logL_projected": float(ll_opt),
    "seed": SEED,
    "projected_fit": {"nu": float(popt[0]), "T_coal": float(popt[1]),
                      "theta": float(theta_proj), "Nref": float(Nref_proj)},
    "adopted_point_fit": {"nu": NU_PUB, "T_coal": T_PUB, "theta": THETA_PUB,
                          "N_ANC": N_ANC_PUB, "N0": N0_PUB, "T_generations": T_YEARS_PUB},
    "gim_sd_linear": {"nu": float(uncerts_lin[0]), "T_coal": float(uncerts_lin[1]),
                      "theta": float(uncerts_lin[2])},
    "gim_sd_log_relative": {"nu": float(uncerts_log[0]), "T_coal": float(uncerts_log[1]),
                            "theta": float(uncerts_log[2])},
    "cov_log_nu_T_theta": [[float(x) for x in row] for row in cov_log],
    "params": {k: summ(derived[k], point_vals[k]) for k in derived},
    "nt_envelope": {"t": [float(x) for x in tgrid], "lo": [float(x) for x in env_lo],
                    "hi": [float(x) for x in env_hi], "median": [float(x) for x in env_med]},
    "notes": ("Blocks = whole chromosomes (LD-respecting). Per-replicate ANGSD/realSFS GL "
              "EM is infeasible at this scale, so the per-chromosome folded SFS point "
              "estimates are resampled and re-summed. Spectra projected to n=%d for fast, "
              "stable moments integration (nu, T, theta are sample-size independent). Point "
              "= adopted full-data growth fit; GIM log-covariance propagated to demographic "
              "units by Monte Carlo centred on that point." % PROJ),
}
out = f"{D}/moments_godambe_uncert.json"
json.dump(results, open(out, "w"), indent=2)
print(f"\nwrote {out}")
for k in ["N_ANC", "N0_present", "T_generations", "expansion_factor"]:
    s = results["params"][k]
    print(f"  {k:16s} point={s['point']:>14,.0f}  95% CI=[{s['ci95'][0]:>13,.0f}, {s['ci95'][1]:>13,.0f}]  CV={s['cv']:.3f}")
