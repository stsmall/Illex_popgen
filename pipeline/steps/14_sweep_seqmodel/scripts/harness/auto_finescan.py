#!/usr/bin/env python3
"""Fine-scale localization + Tajima's-D-sign triage of the 11 BGS-robust autosomal sweep candidates,
using ANGSD windowed_diversity.tsv (10kb: chrom WinCenter pi thetaW Tajima nSites). For each: locate the
pi-min valley bottom in +/-0.5Mb, its D (sign triage), whether the bottom is IN the 100kb call or on the
shoulder, chrom pi-percentile, and the gene at the bottom."""
import numpy as np, collections
WD="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/windowed_diversity.tsv"
GFF="/sietch_colab/data_share/illex/andy_links/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.gff3"
# 11 robust candidates: chrom, call_start, call_end, class, §8e-label
C=[("1",20300001,20400000,"hard","Aplnr"),("1",35900001,36000000,"hard","desert"),
   ("1",85500001,85600000,"hard","eif3h/babam1"),("5",38000001,38100000,"hard","cept1"),
   ("11",10800001,10900000,"hard","desert"),("16",35300001,35400000,"hard","novel/Euprymna"),
   ("25",65000001,65100000,"soft","tRNA-array"),("32",9300001,9400000,"soft","MACROD2"),
   ("32",24100001,24300000,"hard","desert"),("33",9000001,9100000,"hard","Tinagl1"),
   ("35",39000001,39100000,"hard","ZEB2")]
# load windowed diversity per chrom
W=collections.defaultdict(list)
with open(WD) as fh:
    next(fh)
    for ln in fh:
        f=ln.split("\t"); W[f[0]].append((int(f[1]),float(f[2]),float(f[4]),int(f[5])))  # center,pi,D,nSites
for c in W: W[c]=sorted(W[c])
# genes per chrom
G=collections.defaultdict(list)
for ln in open(GFF):
    if ln[0]=="#": continue
    f=ln.split("\t")
    if f[2]!="gene": continue
    nm=f[8].split("Name=")[1].split(";")[0].strip(' "') if "Name=" in f[8] else f[8].split("ID=")[1].split(";")[0]
    G[f[0]].append((int(f[3]),int(f[4]),nm))
def gene_at(c,a,b):
    hits=[nm for s,e,nm in G.get(c,[]) if s<=b and e>=a]
    return hits[0][:42] if hits else "(gene desert)"
print(f"{'cand(§8e)':<22}{'call':<16}{'valley(pi-min)':<15}{'pi':>8}{'pct':>5}{'D':>7}{'pos':>9}  gene@valley  VERDICT")
for c,cs,ce,cl,lab in C:
    arr=W.get(c,[])
    if not arr: print(f"{c}:{lab} — no windows"); continue
    cen=(cs+ce)//2
    reg=[(p,pi,d,ns) for p,pi,d,ns in arr if abs(p-cen)<=500000]
    allpi=sorted(pi for _,pi,_,_ in arr)
    p,pi,d,ns=min(reg,key=lambda x:x[1])   # pi-min window
    pct=int(100*np.searchsorted(allpi,pi)/len(allpi))
    inwin = "IN-call" if cs<=p<=ce else f"shoulder({(p-cen)//1000:+d}kb)"
    g=gene_at(c,p-7500,p+7500)
    # D-sign triage
    if d>=-0.5 and pct<=15: v="ARTIFACT? (deep pi, D>=-0.5)"
    elif d<=-1.5 and pct<=20: v="SWEEP-like (deep pi, D neg)"
    elif pct>40: v="weak (pi not reduced at valley)"
    else: v="intermediate"
    print(f"{c+':'+str(cen//1000000)+'.'+str((cen//100000)%10)+'M '+lab:<22}{str(cs//1000)+'-'+str(ce//1000)+'k':<16}{str(p//1000)+'k':<15}{pi:>8.5f}{pct:>4}%{d:>7.2f}{inwin:>9}  {g}  {v}")
