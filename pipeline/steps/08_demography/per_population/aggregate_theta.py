#!/usr/bin/env python
"""Aggregate per-(division x chromosome) ANGSD thetaStat outputs into per-population
diversity tables. Route A (ANGSD GL, folded SAF, accessible sites) -- the correct
estimator for this low-coverage unphased data; reproduces the genome-wide pi~0.0093,
Tajima's D~-2.

Inputs  : thetas/<div>.<chrom>.win50k.pestPG (50kb non-overlapping windows)
          thetas/<div>.<chrom>.avg.pestPG     (whole-chromosome point estimate)
Outputs (tables/):
  perpop_windows.tsv        long: division,chrom,WinCenter,pi,thetaW,tajd,nSites  (box-plot source)
  perpop_diversity.tsv      one row/division: pi & D point estimate + window medians + N,lat
  perpop_per_chrom.tsv      division x chrom point estimates (supp source)
"""
import glob, os, numpy as np, pandas as pd

D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/per_population"
TH=f"{D}/thetas"; OUT=f"{D}/tables"; os.makedirs(OUT,exist_ok=True)
KARYO="/sietch_colab/data_share/illex/popgen_data/analysis/steps/04_hwe_cline/karyo_by_division.tsv"
cols=["paren","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites"]
MINSITES=5000   # require >=5kb data in a 50kb window

kary=pd.read_csv(KARYO,sep="\t")[["division","N","mean_lat"]].set_index("division")

# ---- windowed (box-plot source) ----
wframes=[]
for f in sorted(glob.glob(f"{TH}/*.win50k.pestPG")):
    b=os.path.basename(f).split(".win50k")[0]; div,chrom=b.split(".",1)
    df=pd.read_csv(f,sep="\t",header=0,names=cols)
    df=df[df["nSites"]>=MINSITES].copy()
    if not len(df): continue
    df["division"]=div; df["chrom"]=chrom
    df["pi"]=df["tP"]/df["nSites"]; df["thetaW"]=df["tW"]/df["nSites"]; df["tajd"]=df["Tajima"]
    wframes.append(df[["division","chrom","WinCenter","pi","thetaW","tajd","nSites"]])
W=pd.concat(wframes,ignore_index=True)
W.to_csv(f"{OUT}/perpop_windows.tsv",sep="\t",index=False)

# ---- whole-chrom point estimates (per div x chrom) ----
aframes=[]
for f in sorted(glob.glob(f"{TH}/*.avg.pestPG")):
    b=os.path.basename(f).split(".avg")[0]; div,chrom=b.split(".",1)
    df=pd.read_csv(f,sep="\t",header=0,names=cols)
    r=df.iloc[0]
    aframes.append({"division":div,"chrom":chrom,"tP":r["tP"],"tW":r["tW"],
                    "Tajima":r["Tajima"],"nSites":r["nSites"],
                    "pi":r["tP"]/r["nSites"],"thetaW":r["tW"]/r["nSites"]})
A=pd.DataFrame(aframes)
A[["division","chrom","pi","thetaW","Tajima","nSites"]].to_csv(f"{OUT}/perpop_per_chrom.tsv",sep="\t",index=False)

# ---- per-population genome(subset) point estimate ----
rows=[]
CHROMS=sorted(A["chrom"].unique(), key=lambda x:int(x) if x.isdigit() else 99)
for div,g in A.groupby("division"):
    pi=g["tP"].sum()/g["nSites"].sum()                       # callable-weighted per-base pi
    tw=g["tW"].sum()/g["nSites"].sum()
    dwt=(g["Tajima"]*g["nSites"]).sum()/g["nSites"].sum()    # nSites-weighted Tajima's D
    wsub=W[W["division"]==div]
    rows.append({"division":div,
        "N_total":int(kary.loc[div,"N"]) if div in kary.index else np.nan,
        "N_used":25 if (div in kary.index and kary.loc[div,"N"]>25) else (int(kary.loc[div,"N"]) if div in kary.index else np.nan),
        "mean_lat":kary.loc[div,"mean_lat"] if div in kary.index else np.nan,
        "pi_point":pi,"thetaW_point":tw,"tajd_point":dwt,
        "pi_median_win":wsub["pi"].median(),"tajd_median_win":wsub["tajd"].median(),
        "n_windows":int(len(wsub)),"n_chroms":int(g["chrom"].nunique())})
S=pd.DataFrame(rows).sort_values("mean_lat").reset_index(drop=True)
S.to_csv(f"{OUT}/perpop_diversity.tsv",sep="\t",index=False)

with open(f"{OUT}/perpop_diversity.provenance.txt","w") as fh:
    fh.write("Per-population diversity (pi, thetaW, Tajima's D) -- ROUTE A (ANGSD GL)\n")
    fh.write("angsd -doSaf 1 -GL 1 -fold (folded SAF) on accessible sites "
             "(analysis/steps/_shared/accessible_sites.sites), realSFS -> saf2theta -> "
             "thetaStat do_stat (whole-chrom + 50kb non-overlapping windows).\n")
    fh.write("Per-division bamlists from baker dedup BAMs; large divisions SUBSAMPLED "
             "to N=25 individuals (seed 42) to bound I/O -- ample for pi/folded-SFS D.\n")
    fh.write(f"Chromosomes: {','.join(CHROMS)} (representative large autosomes; "
             "chr2 inversion, chr42 sex, chrZ excluded).\n")
    fh.write("pi_point/thetaW_point = callable-weighted per-base (sum tP / sum nSites); "
             "tajd_point = nSites-weighted Tajima's D. *_median_win = median over 50kb windows "
             "(>=5kb data). Division 3N (N=1) dropped.\n")

print("windows:",len(W)," per-div-chrom rows:",len(A))
print(S.to_string(index=False))
