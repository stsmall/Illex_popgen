#!/usr/bin/env python
"""Build genome-wide pairwise Hudson FST matrix (per-population) from the existing
pg_gpu per-chromosome results (pggpu_illex/popstats/results/chrom.*.json).

FST is a variance RATIO; it is far more robust to the variants-only / low-coverage
ascertainment bias than absolute pi/D (both numerator and denominator are affected
alike), so route B (pg_gpu, called genotypes) is appropriate for the FST matrix.

Aggregation across chromosomes: callable-bp-weighted mean of per-chromosome Hudson
FST (matches pggpu_illex/popstats/scripts/aggregate.py). Chromosomes used:
autosomes 1-45 EXCLUDING chr2 (segregating inversion -> spurious structure) and
chr42 (sex chromosome). chrZ excluded (different Ne / sex-linked). Division 3N
(N=1) dropped. Written to tables/fst_matrix.tsv (+ long form).
"""
import glob, json, math, os
import numpy as np, pandas as pd

RES="/sietch_colab/data_share/illex/popgen_data/pggpu_illex/popstats/results"
OUT="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/per_population/tables"
os.makedirs(OUT, exist_ok=True)
EXCLUDE_CHROMS={"2","42","Z"}   # inversion, sex chrom, Z
DROP_POPS={"3N"}

data={}
for f in glob.glob(f"{RES}/chrom.*.json"):
    d=json.load(open(f)); c=str(d["chrom"])
    if c in EXCLUDE_CHROMS: continue
    data[c]=d
chroms=sorted(data)
pops=sorted({p for c in chroms for p in data[c]["pops"] if p not in DROP_POPS})
idx={p:i for i,p in enumerate(pops)}
num=np.zeros((len(pops),len(pops))); den=np.zeros((len(pops),len(pops)))
long=[]
for c in chroms:
    call=data[c]["callable"]
    for pair,v in data[c]["pairs"].items():
        a,b=pair.split("__")
        if a in DROP_POPS or b in DROP_POPS: continue
        if a not in idx or b not in idx: continue
        val=v["fst"]
        if not math.isfinite(val): continue
        i,j=idx[a],idx[b]
        num[i,j]+=val*call; num[j,i]+=val*call
        den[i,j]+=call;     den[j,i]+=call
        long.append({"chrom":c,"pop1":a,"pop2":b,"fst":val,"callable":call})
with np.errstate(invalid="ignore"):
    mat=np.where(den>0,num/den,np.nan)
np.fill_diagonal(mat,0.0)
M=pd.DataFrame(mat,index=pops,columns=pops)
M.to_csv(f"{OUT}/fst_matrix.tsv",sep="\t")
pd.DataFrame(long).to_csv(f"{OUT}/fst_per_chrom_long.tsv",sep="\t",index=False)

# provenance
with open(f"{OUT}/fst_matrix.provenance.txt","w") as fh:
    fh.write("Pairwise Hudson FST matrix, per NAFO division (genome-wide)\n")
    fh.write("Source: pg_gpu (route B) per-chromosome divergence.fst on baker SNP "
             "callset (pggpu_illex/popstats/results/chrom.*.json), computed by "
             "pggpu_illex/popstats/scripts/compute_chrom.py (HaplotypeMatrix.from_vcf "
             "-> divergence.fst, Hudson estimator).\n")
    fh.write(f"Chromosomes aggregated (callable-bp-weighted mean): {','.join(chroms)}\n")
    fh.write("EXCLUDED: chr2 (segregating inversion), chr42 (sex chromosome), chrZ.\n")
    fh.write("Dropped population 3N (N=1).\n")
    fh.write("FST is a variance ratio -> robust to the variants-only/low-coverage "
             "ascertainment bias that biases absolute pi/D.\n")

vals=M.values[np.triu_indices(len(pops),1)]
print("pops:",pops)
print("n chroms:",len(chroms))
print(f"FST off-diagonal: min={np.nanmin(vals):.5f} max={np.nanmax(vals):.5f} "
      f"mean={np.nanmean(vals):.5f} median={np.nanmedian(vals):.5f}")
print("frac negative:", float(np.mean(vals<0)))
print(M.round(4).to_string())
