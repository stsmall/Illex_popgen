#!/usr/bin/env python
"""Constraint check on the Z/autosome dS ratio (Illex vs Sthenoteuthis). Lower Z dN alone does not mean
stronger constraint, because a lower Z mutation rate lowers dN and dS together, so genes are matched on
omega = dN/dS (constraint per mutation) rather than dN. Estimates of R = mean dS_Z / mean dS_A and
alpha = (3R-2)/(4-3R):
  all             all genes (as in stheno_ds.summary.json)
  omega_matched   each Z gene against autosomal genes in the same omega decile
  low_constraint  genes above the autosomal median omega only (less constrained)
  gc3_matched     each Z gene against autosomal genes in the same GC3 decile (Illex CDS)
  omega_gc3       matched jointly on omega and GC3 quintiles
CIs: 2,000 bootstrap resamples of Z and autosomal genes within each stratum (Z genes keep their stratum counts).
Usage: alpha_constraint_check.py PER_GENE_TSV ILLEX_CDS_FA OUT_JSON"""
import sys, json, numpy as np, pandas as pd
G, CDS, OUT = sys.argv[1:4]
g = pd.read_csv(G, sep="\t", dtype={"chrom": str})
g = g[~g.excluded & g.dS.notna() & g.dN.notna() & (g.dS > 0)].copy()
want = set(g.tx); seq, k = {}, None
for l in open(CDS):
    if l.startswith(">"): k = l[1:].split()[0]; seq[k] = []
    elif k in want: seq[k].append(l.strip().upper())
def gc3(s):
    s = "".join(s); t = s[2::3]; return (t.count("G") + t.count("C")) / max(len(t), 1)
g["gc3"] = g.tx.map(lambda t: gc3(seq.get(t, [])))
g["omega"] = g.dN / g.dS
g["isZ"] = g.chrom == "Z"
alpha = lambda R: (3 * R - 2) / (4 - 3 * R)
rng = np.random.default_rng(1)
def strata(z, a, keys):
    """per-stratum dS arrays for Z and autosomes (strata = keys), strata present in both."""
    kz = z[keys].astype(str).agg("|".join, axis=1) if keys else pd.Series("all", index=z.index)
    ka = a[keys].astype(str).agg("|".join, axis=1) if keys else pd.Series("all", index=a.index)
    common = sorted(set(kz) & set(ka))
    return [(z.dS.to_numpy()[kz.to_numpy() == c], a.dS.to_numpy()[ka.to_numpy() == c]) for c in common]
def R_of(st):
    """mean dS_Z over the Z-count-weighted mean of stratum autosomal dS means (unweighted if one stratum)."""
    nz = np.array([len(zz) for zz, _ in st]); mz = np.concatenate([zz for zz, _ in st]).mean()
    return mz / np.average([aa.mean() for _, aa in st], weights=nz)
def run(name, z, a, keys=None):
    st = strata(z, a, keys or []); est = R_of(st); bs = []
    for _ in range(2000):
        bs.append(R_of([(zz[rng.integers(0, len(zz), len(zz))], aa[rng.integers(0, len(aa), len(aa))]) for zz, aa in st]))
    ci = np.percentile(bs, [2.5, 97.5])
    return {"n_Z": int(sum(len(zz) for zz, _ in st)), "n_A": int(sum(len(aa) for _, aa in st)), "R": float(est), "R_CI": ci.tolist(),
            "alpha": float(alpha(est)), "alpha_CI": [float(alpha(ci[0])), float(alpha(ci[1]))]}
Z, A = g[g.isZ].copy(), g[~g.isZ].copy()
res = {"omega_median_Z": float(Z.omega.median()), "omega_median_A": float(A.omega.median()),
       "gc3_mean_Z": float(Z.gc3.mean()), "gc3_mean_A": float(A.gc3.mean())}
from scipy.stats import mannwhitneyu
res["MW_omega_p"] = float(mannwhitneyu(Z.omega, A.omega).pvalue); res["MW_gc3_p"] = float(mannwhitneyu(Z.gc3, A.gc3).pvalue)
res["all"] = run("all", Z, A)
edges = np.quantile(A.omega, np.linspace(0, 1, 11)); edges[0], edges[-1] = -np.inf, np.inf
for d in (Z, A): d["om_bin"] = np.digitize(d.omega, edges[1:-1])
res["omega_matched"] = run("om", Z, A, ["om_bin"])
med = A.omega.median()
res["low_constraint"] = run("low", Z[Z.omega > med], A[A.omega > med])
e2 = np.quantile(A.gc3, np.linspace(0, 1, 11)); e2[0], e2[-1] = -np.inf, np.inf
for d in (Z, A): d["gc_bin"] = np.digitize(d.gc3, e2[1:-1])
res["gc3_matched"] = run("gc", Z, A, ["gc_bin"])
e5o = np.quantile(A.omega, np.linspace(0, 1, 6)); e5o[0], e5o[-1] = -np.inf, np.inf
e5g = np.quantile(A.gc3, np.linspace(0, 1, 6)); e5g[0], e5g[-1] = -np.inf, np.inf
for d in (Z, A): d["o5"] = np.digitize(d.omega, e5o[1:-1]); d["g5"] = np.digitize(d.gc3, e5g[1:-1])
res["omega_gc3"] = run("og", Z, A, ["o5", "g5"])
json.dump(res, open(OUT, "w"), indent=1)
for k, v in res.items(): print(k, v)
