#!/usr/bin/env python
"""Aggregate per-(division x chromosome) ANGSD thetaStat outputs into an EXPANDED
windowed diversity table for the manuscript diversity supplement figure.

IDENTICAL method to aggregate_theta.py (Route A: ANGSD GL, folded SAF, accessible
sites; 50kb non-overlapping windows; pi=tP/nSites, thetaW=tW/nSites, tajd=Tajima;
require >=5kb data per window). Only difference: covers the original 3 chroms
(1,13,21) PLUS the added chroms. Writes perpop_windows_expanded.tsv; does NOT touch
the original perpop_windows.tsv.
"""
import glob, os, numpy as np, pandas as pd

D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/per_population"
TH=f"{D}/thetas"; OUT=f"{D}/tables"; os.makedirs(OUT,exist_ok=True)
cols=["paren","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites"]
MINSITES=5000

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
# stable ordering: division, numeric chrom, WinCenter
W["_c"]=W["chrom"].apply(lambda x:int(x) if str(x).isdigit() else 999)
W=W.sort_values(["division","_c","WinCenter"]).drop(columns="_c").reset_index(drop=True)
W.to_csv(f"{OUT}/perpop_windows_expanded.tsv",sep="\t",index=False)

chroms=sorted(W["chrom"].unique(), key=lambda x:int(x) if str(x).isdigit() else 999)
print("chroms:",",".join(chroms),"| total windows:",len(W),"| divisions:",W["division"].nunique())
print("\nPer-chromosome median (pooled over divisions & windows):")
print(f"{'chrom':>6} {'n_win':>7} {'median_pi':>11} {'median_D':>10}")
for c in chroms:
    g=W[W["chrom"]==c]
    print(f"{c:>6} {len(g):>7} {g['pi'].median():>11.6f} {g['tajd'].median():>10.4f}")
print("\nGENOME (all chroms pooled): median_pi=%.6f  median_D=%.4f"%(W['pi'].median(),W['tajd'].median()))
