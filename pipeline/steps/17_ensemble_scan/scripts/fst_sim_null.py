#!/usr/bin/env python
"""
Neutral simulation null for per-SNP geographic (NAFO-division) Weir-Cockerham FST.

The 10 sampled NAFO divisions show NO population structure (global pairwise
Fst ~= 0, PCA panmictic). We build the null explicitly: simulate a SINGLE
panmictic population under the fitted 2-epoch growth demography (moments), draw
the *empirical* per-division sample sizes, and randomly assign the simulated
individuals to the same 10 "divisions" (no real structure). Per-SNP WC FST is
then computed across the 45 division pairs exactly as for the empirical data
(same >=40 obs/count filter; per-SNP mean across pairs with >=20 informative
pairs), reproducing the drift + finite-sample + missing-data noise floor. The
99.9th percentile / max of this null is the outlier threshold.

The empirical callset is ~44% missing, which inflates per-SNP FST sampling
variance; we reproduce that by masking each simulated genotype with its matched
empirical sample's F_MISS.

Efficient design: per-SNP per-POPULATION sufficient statistics (n, allele
count, het count) are computed once for the 10 divisions; the 45 pairwise WC
FSTs are then cheap (NS,)-vector combinations.
"""
import os, sys, time, json
import numpy as np
import msprime

WD  = "/tmp/claude-1003/-sietch-colab-data-share-illex-popgen-data-mkado-illex/6d7afefb-6654-47e0-96f2-bc56ef699fbb/scratchpad/king"
OUT = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/17_ensemble_scan/results/fst_persnp"
LOG = f"{OUT}/sim_progress.log"
os.makedirs(OUT, exist_ok=True)

def log(msg):
    with open(LOG, "a") as fh:
        fh.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")
    print(msg, flush=True)

open(LOG, "w").close()

# ----- fitted demography (moments 2-epoch exponential growth) -----------------
N0, N_ANC, T = 6808096.0, 547928.0, 769519.0
MU, RHO = 3e-9, 0.0     # recomb irrelevant to the MARGINAL per-SNP FST null:
                        # SNPs are treated marginally, so linkage only correlates
                        # neighbouring SNPs, it does not bias the per-SNP FST
                        # distribution. Setting r=0 and using MANY short
                        # independent replicate genealogies gives the correct
                        # marginal null far more cheaply than one recombining ARG.
SEQLEN_SEG = 5e3
TARGET_SNP = 180000     # collect independent trees until >= this many SNPs
SEED = 20240915
alpha = np.log(N0 / N_ANC) / T

# ----- empirical sample -> division -> missingness map ------------------------
div = {}
with open(f"{WD}/pop2.txt") as fh:
    next(fh)
    for ln in fh:
        iid, pop = ln.split()[:2]; div[iid] = pop
fmiss = {}
with open(f"{WD}/miss.smiss") as fh:
    next(fh)
    for ln in fh:
        f = ln.split(); fmiss[f[0]] = float(f[-1])
samps = [s for s in div if s in fmiss]
DIVS = sorted(set(div[s] for s in samps))
order, LAB, FMISS = [], [], []
for d in DIVS:
    for s in sorted(s for s in samps if div[s] == d):
        order.append(s); LAB.append(d); FMISS.append(fmiss[s])
LAB = np.array(LAB); FMISS = np.array(FMISS); NIND = len(LAB)
div_idx = {d: np.where(LAB == d)[0] for d in DIVS}
PAIRS = [(i, j) for i in range(len(DIVS)) for j in range(i+1, len(DIVS))]
assert len(PAIRS) == 45
log(f"[map] {NIND} indiv, {len(DIVS)} divisions, mean F_MISS={FMISS.mean():.3f}")
log(f"[map] sizes: {dict((d,int((LAB==d).sum())) for d in DIVS)}")
log(f"[demog] alpha={alpha:.4e} check N0*exp(-aT)={N0*np.exp(-alpha*T):.0f}")

# ----- simulate ---------------------------------------------------------------
def sim_block(seed):
    dem = msprime.Demography()
    dem.add_population(name="A", initial_size=N0, growth_rate=alpha)
    dem.add_population_parameters_change(time=T, population="A",
                                         growth_rate=0, initial_size=N_ANC)
    ts = msprime.sim_ancestry(samples=NIND, demography=dem,
                              sequence_length=SEQLEN_SEG, recombination_rate=RHO,
                              ploidy=2, random_seed=seed)
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=seed)
    H = ts.genotype_matrix()
    return (H[:, 0::2] + H[:, 1::2]).astype(np.int8)

t0 = time.time(); rng = np.random.default_rng(SEED)
blocks = []; nsofar = 0; k = 0
while nsofar < TARGET_SNP:
    b = sim_block(SEED + 1 + k)
    blocks.append(b); nsofar += b.shape[0]; k += 1
    if k % 100 == 0:
        log(f"[sim] {k} indep trees, {nsofar} SNPs ({time.time()-t0:.1f}s)")
G = np.concatenate(blocks, axis=0); del blocks
NS = G.shape[0]; N_SEG = k
log(f"[sim] {NS} SNPs from {k} independent genealogies x {NIND} indiv in {time.time()-t0:.1f}s")

# ----- inject empirical per-sample missingness (code missing as -1) -----------
for i in range(NIND):
    if FMISS[i] > 0:
        G[rng.random(NS) < FMISS[i], i] = -1
log(f"[miss] injected missing fraction = {np.mean(G < 0):.3f}")

# ----- per-SNP per-POP sufficient stats (computed once) -----------------------
# transpose so each individual's genotypes are contiguous -> fast row groups
GT = np.ascontiguousarray(G.T)            # (NIND, NS) int8
Npop  = np.zeros((len(DIVS), NS), np.float64)
ACpop = np.zeros((len(DIVS), NS), np.float64)
HETpop= np.zeros((len(DIVS), NS), np.float64)
t1 = time.time()
for pi, d in enumerate(DIVS):
    sub = GT[div_idx[d]]                   # (n_p, NS) contiguous rows
    miss = sub < 0
    Npop[pi]  = (~miss).sum(0)
    ACpop[pi] = np.where(miss, 0, sub).sum(0)
    HETpop[pi]= (sub == 1).sum(0)
    log(f"[stats] pop {d} ({pi+1}/10) done ({time.time()-t1:.1f}s)")
del GT, G

# ----- 45 pairwise WC FST (cheap vector combinations) -------------------------
pooled = []
psum = np.zeros(NS); pcnt = np.zeros(NS); pmax = np.zeros(NS)
t2 = time.time()
for k, (i, j) in enumerate(PAIRS):
    n1 = Npop[i]; n2 = Npop[j]
    ac1 = ACpop[i]; ac2 = ACpop[j]
    h1c = HETpop[i]; h2c = HETpop[j]
    obs = n1 + n2
    with np.errstate(divide='ignore', invalid='ignore'):
        p1 = ac1/(2*n1); p2 = ac2/(2*n2)
        hb1 = h1c/n1;    hb2 = h2c/n2
        r = 2.0; N = n1 + n2; nbar = N/r
        nc = (N - (n1*n1 + n2*n2)/N) / (r-1)
        pbar = (n1*p1 + n2*p2)/N
        hbar = (n1*hb1 + n2*hb2)/N
        s2 = (n1*(p1-pbar)**2 + n2*(p2-pbar)**2) / ((r-1)*nbar)
        a = (nbar/nc)*(s2 - (1/(nbar-1))*(pbar*(1-pbar) - ((r-1)/r)*s2 - hbar/4))
        b = (nbar/(nbar-1))*(pbar*(1-pbar) - ((r-1)/r)*s2 - ((2*nbar-1)/(4*nbar))*hbar)
        c = hbar/2
        den = a + b + c
        fst = np.where(den != 0, a/den, np.nan)
    bad = (n1 < 2) | (n2 < 2) | ~np.isfinite(fst)
    fst = np.where(bad, np.nan, fst)
    keep = (obs >= 40) & np.isfinite(fst)
    # FST is defined on [0,1]; the WC *estimator* can stray outside under small
    # effective n (rare variants in the smallest divisions x 44% missing) -> clip
    # to the biological range, as is standard for per-SNP FST (applied identically
    # to the empirical data downstream).
    fclip = np.clip(fst, 0, 1.0)
    pooled.append(fclip[keep].astype(np.float32))
    psum[keep] += fclip[keep]; pcnt[keep] += 1
    pmax = np.where(keep & (fclip > pmax), fclip, pmax)
    if k % 9 == 8:
        log(f"[fst] pair {k+1}/45 done ({time.time()-t2:.1f}s)")
pooled = np.concatenate(pooled)
log(f"[fst] {len(pooled):,} obs-filtered per-SNP-per-pair values ({time.time()-t2:.1f}s)")

ok = pcnt >= 20
meanf = psum[ok]/pcnt[ok]
def q(a,p): return float(np.quantile(a,p))
summary = {
  "n_sim_snp_perpair": int(len(pooled)), "n_sim_snp_persnp": int(ok.sum()),
  "n_sim_segments": N_SEG, "n_indiv": int(NIND),
  "demography": {"N0":N0,"N_ANC":N_ANC,"T_gen":T,"mu":MU,"recomb":RHO,"growth_rate":alpha},
  "pooled_perpair": {"median":float(np.median(pooled)),"mean":float(pooled.mean()),
     "q99":q(pooled,.99),"q999":q(pooled,.999),
     "frac_gt_0.25":float(np.mean(pooled>0.25)),"frac_gt_0.5":float(np.mean(pooled>0.5)),
     "max":float(pooled.max())},
  "persnp_mean": {"median":float(np.median(meanf)),"q99":q(meanf,.99),
     "q999":q(meanf,.999),"q9999":q(meanf,.9999),"max":float(meanf.max())},
}
log("[SUMMARY] " + json.dumps(summary, indent=2))
np.save(f"{OUT}/sim_pooled_fst.npy", pooled)
np.save(f"{OUT}/sim_persnp_meanfst.npy", meanf.astype(np.float32))
with open(f"{OUT}/sim_null_summary.json","w") as fh: json.dump(summary, fh, indent=2)
log("[done] wrote sim outputs")
