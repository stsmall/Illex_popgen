"""Permutation null for inversion-gene GO enrichment: draw K random gene sets of the
same size (29 GO-annotated genes) and, at each QQ rank, take the 2.5-97.5% band of the
sorted -log10(p) vectors. This captures the GO-hierarchy + small-count inflation, so the
band is the RIGHT null. Also counts nominal-significant terms per draw vs observed.
Writes perm_null_band.tsv (rank, exp, lo, hi) + perm_null_counts.txt."""
import numpy as np, pandas as pd
from scipy.stats import hypergeom
from scipy import sparse
A = "/sietch_colab/data_share/illex/popgen_data"
GENEGO = f"{A}/analysis/steps/17_ensemble_scan/results/fst_persnp/gene_go.tsv"
O = f"{A}/analysis/steps/03_karyotype/inversion_content"
K = 500
rng = np.random.default_rng(0)

gg = pd.read_csv(GENEGO, sep="\t", header=None, names=["gene", "go"])
genes = gg.gene.unique(); terms = gg.go.unique()
gidx = {g: i for i, g in enumerate(genes)}; tidx = {t: i for i, t in enumerate(terms)}
r = gg.gene.map(gidx).values; c = gg.go.map(tidx).values
M = sparse.csr_matrix((np.ones(len(gg)), (r, c)), shape=(len(genes), len(terms)))  # gene x term
Ng = len(genes); nt = np.asarray(M.sum(0)).ravel()  # genes per term
N_INV = 29

# observed p-values already computed -> load, align to terms present here
obs = pd.read_csv(f"{O}/inversion_GO_enrichment_allcat.tsv", sep="\t")
obs_y = np.sort(-np.log10(obs.p.clip(lower=1e-300).values))[::-1]
n_obs_nominal = int((obs.p < 0.05).sum())


def enrich_neglogp(sel_rows):
    a = np.asarray(M[sel_rows].sum(0)).ravel().astype(int)   # overlap per term
    keep = a >= 1
    p = np.ones(len(terms))
    p[keep] = hypergeom.sf(a[keep] - 1, Ng, nt[keep].astype(int), N_INV)
    p = p[a >= 2]                                            # match observed's a>=2 rule
    return np.sort(-np.log10(np.clip(p, 1e-300, 1)))[::-1], int((p < 0.05).sum())


null_sorted, null_counts = [], []
for k in range(K):
    sel = rng.choice(Ng, N_INV, replace=False)
    y, nc = enrich_neglogp(sel)
    null_sorted.append(y); null_counts.append(nc)

# align lengths (# terms with a>=2 varies) -> pad to common max length with 0
L = max(len(v) for v in null_sorted + [obs_y])
def pad(v): return np.concatenate([v, np.zeros(L - len(v))])
NS = np.vstack([pad(v) for v in null_sorted])
exp = NS.mean(0); lo = np.percentile(NS, 2.5, 0); hi = np.percentile(NS, 97.5, 0)
oy = pad(obs_y)
band = pd.DataFrame({"rank": np.arange(1, L + 1), "obs": oy, "exp": exp, "lo": lo, "hi": hi})
band.to_csv(f"{O}/perm_null_band.tsv", sep="\t", index=False)
nc = np.array(null_counts)
with open(f"{O}/perm_null_counts.txt", "w") as fh:
    frac = float((nc >= n_obs_nominal).mean())
    fh.write(f"observed_nominal={n_obs_nominal}\nnull_nominal_mean={nc.mean():.1f}\n"
             f"null_nominal_sd={nc.std():.1f}\nnull_nominal_range={nc.min()}-{nc.max()}\n"
             f"p_observed_ge_null={frac:.3f}\nK={K}\n")
# fraction of observed points outside the permutation band
out_hi = float((oy > hi + 1e-9).mean())
print(f"DONE_PERM observed_nominal={n_obs_nominal} null_mean={nc.mean():.1f}+-{nc.std():.1f} "
      f"(range {nc.min()}-{nc.max()}) p={frac:.3f} | obs points above band: {100*out_hi:.1f}%")
