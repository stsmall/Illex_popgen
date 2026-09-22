"""Between-division Hudson FST (ratio of averages, all 45 division pairs pooled) in 100-kb windows across
four elevated-pi blocks and 5 Mb of flank on each side. Input: per-division AC/AN from bcftools +fill-tags."""
import pandas as pd, numpy as np, itertools, sys
BL={"1":(23.98e6,31.28e6),"30":(17.3e6,18.5e6),"23":(0.1e6,5.6e6),"38":(61.5e6,68.0e6)}
divs=sorted(set(l.split("\t")[1].strip() for l in open("sample_groups.txt")))
out=[]
for c,(s,e) in BL.items():
    cols=["pos"]+[f"{k}_{d}" for d in divs for k in ("AC","AN")]
    t=pd.read_csv(f"af_chr{c}.tsv",sep="\t",header=None,usecols=range(len(cols)),names=cols)
    num=np.zeros(len(t)); den=np.zeros(len(t))
    for a,b in itertools.combinations(divs,2):
        n1,n2=t[f"AN_{a}"].values.astype(float),t[f"AN_{b}"].values.astype(float)
        ok=(n1>=10)&(n2>=10)
        p1=np.where(ok,t[f"AC_{a}"]/np.maximum(n1,1),0); p2=np.where(ok,t[f"AC_{b}"]/np.maximum(n2,1),0)
        num+=np.where(ok,(p1-p2)**2-p1*(1-p1)/np.maximum(n1-1,1)-p2*(1-p2)/np.maximum(n2-1,1),0)
        den+=np.where(ok,p1*(1-p2)+p2*(1-p1),0)
    t["win"]=(t.pos//100000)*100000; t["num"]=num; t["den"]=den
    g=t.groupby("win").agg(num=("num","sum"),den=("den","sum"),n=("pos","size")).reset_index()
    g["fst"]=g.num/g.den; g["chrom"]=c; g["region"]=np.where((g.win>=s)&(g.win<e),"block","flank")
    out.append(g)
    ins=g[g.region=="block"]; fl=g[g.region=="flank"]
    print(f"chr{c} block {s/1e6:.1f}-{e/1e6:.1f}: FST block {ins.num.sum()/ins.den.sum():+.5f} (window median {ins.fst.median():+.5f}, max {ins.fst.max():.4f}) | flank {fl.num.sum()/fl.den.sum():+.5f} (median {fl.fst.median():+.5f}, max {fl.fst.max():.4f}) | SNPs/100kb block {ins.n.median():.0f} flank {fl.n.median():.0f}")
pd.concat(out).to_csv("block_fst_windows.tsv",sep="\t",index=False)
