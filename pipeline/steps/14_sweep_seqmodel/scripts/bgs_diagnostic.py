#!/usr/bin/env python3
"""Genome-wide BGS diagnostic: does windowed diversity track CDS density / recombination?

Uses EXISTING ANGSD win50k pestPG (correct invariant-aware diversity via GLs) — no new
sim. BGS signature: pi NEGATIVELY correlated with CDS density and POSITIVELY with
recombination; Tajima's D more negative where CDS-dense / low-recomb.
"""
import os, subprocess, numpy as np, pandas as pd
A = "/sietch_colab/data_share/illex/popgen_data/analysis"
THETA = f"{A}/steps/08_demography/thetas"
CDS = "/sietch_colab/data_share/illex/popgen_data/degenotate_illex/illex.cds_sites.merge.bed"
RMAP = f"{A}/steps/11_relernn/run_male_auto/proj/male.kept.PREDICT.BSCORRECTED.txt"
BB = "/home/ssmall/miniforge3/envs/bioinfo-buddy/bin"
OUT = f"{A}/steps/14_sweep_seqmodel/results"
os.makedirs(OUT, exist_ok=True)
EXCLUDE = {"2", "Z"}          # chr2 inversion, Z sex chromosome
MIN_NSITES = 5000            # require >=5k callable sites in a 50kb window

# ---- 1. parse pestPG (all autosomes) ----
rows = []
for f in sorted(os.listdir(THETA)):
    if not f.endswith(".win50k.pestPG"):
        continue
    chrom = f.split(".")[1]
    if chrom in EXCLUDE:
        continue
    for line in open(f"{THETA}/{f}"):
        if line.startswith("#") or not line.strip():
            continue
        c = line.rstrip("\n").split("\t")
        # cols: [0]blob [1]Chr [2]WinCenter [3]tW [4]tP [8]Tajima [13]nSites
        center = int(c[2]); tP = float(c[4]); tajd = float(c[8]); nsites = int(c[13])
        if nsites < MIN_NSITES:
            continue
        rows.append((chrom, center - 25000, center + 25000, tP / nsites, tajd, nsites))
df = pd.DataFrame(rows, columns=["chrom", "start", "stop", "pi", "tajD", "nSites"])
df["start"] = df["start"].clip(lower=0)
print(f"windows (nSites>={MIN_NSITES}, excl {sorted(EXCLUDE)}): {len(df)}  over {df.chrom.nunique()} chroms")
print(f"genome-wide pi mean={df.pi.mean():.5f}  tajD mean={df.tajD.mean():.3f}")

# ---- 2. CDS density per window (bedtools coverage) ----
wbed = f"{OUT}/_windows.bed"
df.sort_values(["chrom", "start"]).to_csv(wbed, sep="\t", header=False, index=False,
                                          columns=["chrom", "start", "stop"])
cov = subprocess.run([f"{BB}/bedtools", "coverage", "-a", wbed, "-b", CDS],
                     capture_output=True, text=True, check=True).stdout
# coverage cols: chrom start stop  n_overlap  covered_bp  win_len  frac_covered
covmap = {}
for line in cov.strip().split("\n"):
    p = line.split("\t")
    covmap[(p[0], int(p[1]))] = float(p[-1])   # fraction covered by CDS
df["cds_frac"] = [covmap.get((r.chrom, r.start), np.nan) for r in df.itertuples()]

# ---- 3. recombination per window (ReLERNN male_auto, midpoint lookup) ----
mp = {}
for line in open(RMAP):
    if line.startswith("chrom"):
        continue
    p = line.rstrip("\n").split("\t"); ch = p[0].replace("b", "").replace("'", "")
    mp.setdefault(ch, []).append((int(p[1]), int(p[2]), float(p[4])))
def recomb_for(ch, start, stop):
    ivs = mp.get(ch);  mid = (start + stop) // 2
    if not ivs:
        return np.nan
    for s, e, r in ivs:
        if s <= mid < e:
            return r
    return np.nan
df["recomb"] = [recomb_for(r.chrom, r.start, r.stop) for r in df.itertuples()]

d = df.dropna(subset=["pi", "cds_frac", "recomb"]).copy()
print(f"windows with pi+cds+recomb: {len(d)}")

# ---- 4. correlations ----
def spearman(a, b):
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    return np.corrcoef(ra, rb)[0, 1]
print("\n=== BGS signature tests (genome-wide, all autosomes) ===")
print(f"Spearman(pi,   CDS density) = {spearman(d.pi, d.cds_frac):+.3f}   (BGS: NEGATIVE)")
print(f"Spearman(pi,   recomb)      = {spearman(d.pi, d.recomb):+.3f}   (BGS: POSITIVE)")
print(f"Spearman(tajD, CDS density) = {spearman(d.tajD, d.cds_frac):+.3f}   (BGS: NEGATIVE)")
print(f"Spearman(tajD, recomb)      = {spearman(d.tajD, d.recomb):+.3f}   (BGS: POSITIVE)")
print(f"Spearman(cds_frac, recomb)  = {spearman(d.cds_frac, d.recomb):+.3f}   (collinearity check)")

# standardized multiple regression pi ~ cds_frac + recomb
X = np.column_stack([np.ones(len(d)),
                     (d.cds_frac - d.cds_frac.mean()) / d.cds_frac.std(),
                     (d.recomb - d.recomb.mean()) / d.recomb.std()])
y = (d.pi - d.pi.mean()) / d.pi.std()
beta, *_ = np.linalg.lstsq(X, y.values, rcond=None)
print(f"\nstandardized regression  pi ~ CDS + recomb:  beta_CDS={beta[1]:+.3f}  beta_recomb={beta[2]:+.3f}")

# ---- 5. binned table by CDS-density quartile ----
d["cds_q"] = pd.qcut(d.cds_frac.rank(method="first"), 4, labels=["Q1(low)", "Q2", "Q3", "Q4(high)"])
print("\nCDS-density quartile   mean_cds_frac   mean_pi     mean_tajD    n")
for q, g in d.groupby("cds_q", observed=True):
    print(f"  {str(q):16s}   {g.cds_frac.mean():.4f}      {g.pi.mean():.5f}   {g.tajD.mean():+.3f}   {len(g)}")

d.to_csv(f"{OUT}/bgs_diagnostic_windows.tsv", sep="\t", index=False)
print(f"\nwrote {OUT}/bgs_diagnostic_windows.tsv")
