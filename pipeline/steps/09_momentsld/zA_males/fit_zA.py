#!/usr/bin/env python
"""rho_Z / rho_A for ZZ males from the momentsLD statistics in ld_regions_{A,Z}.pkl, two ways.
(1) momentsLD fits (constant size and exponential growth) to each set. Both use the same flat map
(1e-8 per bp), so the fitted Ne is proportional to rho per bp and Ne_Z / Ne_A = rho_Z / rho_A.
(2) LD decay: sigma_D^2 = DD / pi2 per distance bin, fitted to the Ohta-Kimura expectation
(10 + rho d) / (22 + 13 rho d + (rho d)^2) plus a constant sample-size floor, giving rho per bp directly.
(3) Curve matching: the distance stretch k that maps the Z sigma_D^2 curve onto the autosomal one.
CIs from 200 bootstrap resamples of regions. Writes fit_zA.json."""
import pickle, json, numpy as np, moments.LD
from scipy.optimize import curve_fit
MD = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/09_momentsld/zA_males"
def load(g):
    d = pickle.load(open(f"{MD}/ld_regions_{g}.pkl", "rb")); return d
dA, dZ = load("A"), load("Z")
rb = np.array(dA["r_bins"]); bp = np.array(dA["bp_bins"], float)
mid = (bp[:-1] + bp[1:]) / 2
def momfit(regions):
    mv = moments.LD.Parsing.bootstrap_data(regions); data = [mv["means"], mv["varcovs"]]
    out = {}
    for name, f, p0, lb, ub in (("constant", moments.LD.Demographics1D.snm, [1e5], [1e2], [1e10]),
                                ("growth", moments.LD.Demographics1D.growth, [2, 0.1, 1e5], [1e-3, 1e-4, 1e2], [1e4, 10, 1e10])):
        try:
            opt, LL = moments.LD.Inference.optimize_log_fmin(p0, data, [f], rs=rb, lower_bound=lb, upper_bound=ub, verbose=0, maxiter=500)
            out[name] = {"opt": list(map(float, opt)), "LL": float(LL), "Ne": float(opt[-1]), "at_bound": bool(opt[-1] > 0.99 * ub[-1] or opt[-1] < 1.01 * lb[-1])}
        except Exception as e:
            out[name] = {"error": str(e)}
    return out
def sigma(regions):
    """sigma_D^2 per bin: pooled DD over pooled pi2 across regions (stats order from moments.LD names)."""
    names = regions[next(iter(regions))]["stats"][0]
    iDD = names.index("DD_0_0"); iPi = names.index("pi2_0_0_0_0")
    sums = np.array([[r["sums"][b] for b in range(len(rb) - 1)] for r in regions.values()])  # region x bin x stat
    return sums[:, :, iDD], sums[:, :, iPi]
def okfit(DD, PI):
    y = DD.sum(0) / PI.sum(0)
    f = lambda d, rho, c: (10 + rho * d) / (22 + 13 * rho * d + (rho * d) ** 2) + c
    p, _ = curve_fit(f, mid, y, p0=[0.01, 0.003], bounds=([1e-6, 0], [10, 1]), maxfev=20000)
    return p, y
res = {"bp_bins": bp.tolist(), "n_regions": {"A": len(dA["regions"]), "Z": len(dZ["regions"])}}
res["momentsLD"] = {"A": momfit(dA["regions"]), "Z": momfit(dZ["regions"])}
for m in ("constant", "growth"):
    a, z = res["momentsLD"]["A"].get(m, {}), res["momentsLD"]["Z"].get(m, {})
    if "Ne" in a and "Ne" in z: res["momentsLD"][f"ratio_{m}"] = z["Ne"] / a["Ne"]
DDa, PIa = sigma(dA["regions"]); DDz, PIz = sigma(dZ["regions"])
(pa, ya), (pz, yz) = okfit(DDa, PIa), okfit(DDz, PIz)
rng = np.random.default_rng(1); br = []
for _ in range(200):
    ia = rng.integers(0, len(DDa), len(DDa)); iz = rng.integers(0, len(DDz), len(DDz))
    try: br.append(okfit(DDz[iz], PIz[iz])[0][0] / okfit(DDa[ia], PIa[ia])[0][0])
    except Exception: pass
res["ld_decay"] = {"rho_per_bp_A": pa[0], "floor_A": pa[1], "rho_per_bp_Z": pz[0], "floor_Z": pz[1],
                   "ratio": pz[0] / pa[0], "ratio_CI": np.percentile(br, [2.5, 97.5]).tolist(), "n_boot": len(br),
                   "sigmaD2_A": ya.tolist(), "sigmaD2_Z": yz.tolist()}
# (3) curve matching: the k that best maps sigma_D^2_Z(d) onto sigma_D^2_A(k d) (log-log interpolation).
# If both curves are functions of rho*d with the same sample size and history shape, k = rho_Z / rho_A.
def kmatch(yz, ya):
    lx = np.log(mid); ks = np.logspace(-1, 1, 801); best = (np.inf, np.nan)
    for k in ks:
        x = lx + np.log(k); ok = (x >= lx[0]) & (x <= lx[-1])
        if ok.sum() < 4: continue
        e = np.mean((np.log(yz[ok]) - np.interp(x[ok], lx, np.log(ya))) ** 2)
        if e < best[0]: best = (e, k)
    return best[1]
k0 = kmatch(yz, ya); bk = []
for _ in range(200):
    ia = rng.integers(0, len(DDa), len(DDa)); iz = rng.integers(0, len(DDz), len(DDz))
    bk.append(kmatch(DDz[iz].sum(0) / PIz[iz].sum(0), DDa[ia].sum(0) / PIa[ia].sum(0)))
res["curve_match"] = {"ratio": k0, "ratio_CI": np.nanpercentile(bk, [2.5, 97.5]).tolist()}
json.dump(res, open(f"{MD}/fit_zA.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
