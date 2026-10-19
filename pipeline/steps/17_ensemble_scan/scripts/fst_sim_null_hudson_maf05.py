#!/usr/bin/env python
"""Neutral null for per-SNP geographic Hudson FST (MAF > 0.05), matching hudson_maf05_agg.py.

A single panmictic population under the fitted moments growth demography is simulated
with msprime. Simulated diploids are assigned to the 10 NAFO divisions at the empirical
sample sizes (3N, n = 1, excluded), each simulated individual takes the genotype
missingness of its matched empirical individual (pooled over chromosomes from the
plink2 --missing output of run_gw_maf05.sh), and SNPs are kept at pooled MAF > 0.05.
Per pair: Hudson numerator and denominator (Bhatia et al. 2013) with the same pair
filters as the empirical data; per SNP: ratio of averages over >= 20 informative pairs.
Independent short non-recombining segments give the correct marginal per-SNP null.
Saves per-SNP null FST and pooled MAF (for frequency-matched comparison).
Usage: fst_sim_null_hudson_maf05.py OUT_DIR TARGET_SNPS N_WORKERS
"""
import sys, os, glob, json, time, itertools
import numpy as np, msprime
from multiprocessing import Pool

OUT, TARGET, NW = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
HERE = os.path.dirname(os.path.abspath(__file__))
N0, N_ANC, T, MU = 6808096.0, 547928.0, 769519.0, 3e-9
ALPHA = np.log(N0 / N_ANC) / T
SEG, PER_BATCH, SEED = 5e3, 40, 20240915

div = {}
for ln in open(f"{HERE}/../config/divisions_632.tsv"):
    if ln.startswith("#"): continue
    iid, pop = ln.split()[:2]; div[iid] = pop
miss_ct, obs_ct = {}, {}
for f in glob.glob(f"{OUT}/perchrom/miss_*.smiss"):
    for ln in open(f):
        if ln.startswith("#"): continue
        x = ln.split(); miss_ct[x[0]] = miss_ct.get(x[0], 0) + int(x[1]); obs_ct[x[0]] = obs_ct.get(x[0], 0) + int(x[2])
IIDS = [i for i in div if i in obs_ct]
DIVS = sorted(set(div[i] for i in IIDS))
LAB = np.array([div[i] for i in IIDS]); FMISS = np.array([miss_ct[i] / obs_ct[i] for i in IIDS])
NIND = len(IIDS); IDX = [np.where(LAB == d)[0] for d in DIVS]; PAIRS = list(itertools.combinations(range(len(DIVS)), 2))
assert len(PAIRS) == 45, len(PAIRS)

def batch(k):
    rng = np.random.default_rng(SEED + 7919 * k)
    dem = msprime.Demography(); dem.add_population(name="A", initial_size=N0, growth_rate=ALPHA)
    dem.add_population_parameters_change(time=T, population="A", growth_rate=0, initial_size=N_ANC)
    gs = []
    for s in range(PER_BATCH):
        ts = msprime.sim_ancestry(samples=NIND, demography=dem, sequence_length=SEG, ploidy=2,
                                  random_seed=SEED + 100000 * k + s + 1)
        ts = msprime.sim_mutations(ts, rate=MU, random_seed=SEED + 100000 * k + s + 1)
        H = ts.genotype_matrix()
        if H.shape[0]: gs.append((H[:, 0::2] + H[:, 1::2]).astype(np.int8))
    if not gs: return np.zeros(0, np.float32), np.zeros(0, np.float32), 0
    G = np.concatenate(gs).T.copy()                       # (NIND, NS)
    G[rng.random(G.shape) < FMISS[:, None]] = -1           # empirical missingness
    obs = (G >= 0); A = np.where(obs, G, 0)
    a_all = A.sum(0); n_all = 2 * obs.sum(0)
    p_all = np.where(n_all > 0, a_all / np.maximum(n_all, 1), 0)
    maf = np.minimum(p_all, 1 - p_all)
    m = maf > 0.05                                         # pooled MAF > 0.05
    maf = maf[m]
    A, obs = A[:, m], obs[:, m]
    Ac = np.vstack([A[ix].sum(0) for ix in IDX]).astype(np.float64)
    Nc = np.vstack([2 * obs[ix].sum(0) for ix in IDX]).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        P = np.where(Nc > 0, Ac / Nc, np.nan)
        ns = np.zeros(Ac.shape[1]); ds = np.zeros(Ac.shape[1]); npair = np.zeros(Ac.shape[1], int)
        for i, j in PAIRS:
            p1, p2, n1, n2 = P[i], P[j], Nc[i], Nc[j]
            num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
            den = p1 * (1 - p2) + p2 * (1 - p1)
            ok = (n1 / 2 + n2 / 2 >= 40) & (n1 >= 2) & (n2 >= 2) & np.isfinite(num) & (den > 0)
            ns += np.where(ok, num, 0); ds += np.where(ok, den, 0); npair += ok
        keep = npair >= 20
    return (ns[keep] / ds[keep]).astype(np.float32), maf[keep].astype(np.float32), int(G.shape[1])

if __name__ == "__main__":
    os.makedirs(f"{OUT}/null", exist_ok=True)
    t0 = time.time(); vals = []; mafs = []; n = 0; nraw = 0; k = 0
    with Pool(NW) as pool:
        while n < TARGET:
            for v, mf, r in pool.imap_unordered(batch, range(k, k + NW * 4)):
                vals.append(v); mafs.append(mf); n += len(v); nraw += r
            k += NW * 4
            print(f"[{time.time()-t0:7.0f}s] batches {k}  raw SNPs {nraw:,}  kept {n:,}", flush=True)
    f = np.concatenate(vals)
    np.save(f"{OUT}/null/null_hudson_maf05.npy", f)
    np.save(f"{OUT}/null/null_maf.npy", np.concatenate(mafs))
    q = lambda p: float(np.quantile(f, p))
    s = {"n_null_snps": int(len(f)), "n_raw_sim_snps": int(nraw), "n_indiv": NIND, "mean_fmiss": float(FMISS.mean()),
         "median": q(.5), "q99": q(.99), "q999": q(.999), "q9999": q(.9999), "q99999": q(.99999), "max": float(f.max()),
         "demography": {"N0": N0, "N_ANC": N_ANC, "T": T, "mu": MU}}
    json.dump(s, open(f"{OUT}/null/null_summary.json", "w"), indent=1)
    print(json.dumps(s, indent=1))
