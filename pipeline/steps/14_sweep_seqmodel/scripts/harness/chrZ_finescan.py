#!/usr/bin/env python3
"""Genome-wide chrZ 10kb pi/Tajima-D fine scan (illex 330 males), then localize each chrZ diploSHIC
candidate: valley bottom (pi-min), depth, D, chrZ pi-percentile, and overlapping gene(s)."""
import sys, numpy as np, allel
VCF=sys.argv[1]; GFF="/sietch_colab/data_share/illex/andy_links/Illex_F24.gene_lnc_pseudo.func.fix.sq3.FINAL.v2.gff3"
W=10000; ZLEN=38234284
print("loading chrZ males...", file=sys.stderr)
c=allel.read_vcf(VCF, fields=['variants/POS','calldata/GT'])
pos=c['variants/POS']; gt=allel.GenotypeArray(c['calldata/GT']); ac=gt.count_alleles()
pi,win,_,_=allel.windowed_diversity(pos,ac,size=W,start=1,stop=ZLEN)
d,_,_=allel.windowed_tajima_d(pos,ac,size=W,start=1,stop=ZLEN)
wc=(win[:,0]+win[:,1])/2
finite=np.isfinite(pi)
pi_sorted=np.sort(pi[finite])
def pctile(v):
    return int(100*np.searchsorted(pi_sorted,v)/len(pi_sorted))
# genes with names
genes=[]
for line in open(GFF):
    if line[0]=="#": continue
    f=line.rstrip("\n").split("\t")
    if f[0]!="Z" or f[2]!="gene": continue
    nm=f[8];
    if "Name=" in nm: nm=nm.split("Name=")[1].split(";")[0].strip(' "')
    else: nm=nm.split("ID=")[1].split(";")[0]
    genes.append((int(f[3]),int(f[4]),nm))
def genes_in(a,b):
    return [f"{nm}({s//1000}-{e//1000}kb)" for s,e,nm in genes if s<=b and e>=a][:4]

cands=[(12750000,"hard",0.95),(23650000,"hard",0.97),(24550000,"soft",0.92),
       (25450000,"soft",0.46),(30750000,"soft",0.78),(31450000,"hard",0.63)]
print(f"chrZ 10kb pi: median={np.median(pi_sorted):.5f}  (n={len(pi_sorted)} win)\n")
print(f"{'candidate':<12}{'class/P':<11}{'valley(pi-min)':<18}{'pi_min':>8}{'pctile':>7}{'D':>7}{'depth':>7}  target gene(s)")
FL=500000
for center,cl,p in cands:
    sel=(wc>=center-FL)&(wc<=center+FL)&finite
    idx=np.where(sel)[0]
    mi=idx[np.argmin(pi[idx])]
    flank=np.median(pi[sel])
    depth=pi[mi]/flank if flank else np.nan
    tg=genes_in(win[mi,0]-15000,win[mi,1]+15000)
    vb=f"Z:{win[mi,0]/1e6:.2f}-{win[mi,1]/1e6:.2f}"
    print(f"Z:{center/1e6:<9.2f}{cl+'/'+str(p):<11}{vb:<18}{pi[mi]:>8.5f}{pctile(pi[mi]):>6}%{d[mi]:>7.2f}{depth:>7.2f}  {', '.join(tg) if tg else '(gene desert)'}")
