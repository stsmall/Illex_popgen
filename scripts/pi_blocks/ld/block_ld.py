"""LD test for cryptic inversions in the elevated-pi blocks. An inversion segregating at appreciable frequency
produces long-range LD spanning the whole block (chr2 inversion: positive control). A fast-evolving,
relaxed-constraint region produces ordinary LD decay. For each region (block +/- 1/3 of its length of flank):
common SNPs (MAF>=0.1, <=30% missing; callset per-site missingness is ~40%), thinned to <=1500, genotype r^2 (0/1/2 dosage, pairwise complete).
Outputs r^2 matrix per region and a summary of mean r^2 by distance class inside vs outside the block."""
import sys, subprocess, numpy as np, pandas as pd, io
BCF="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"
V="/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/{c}/variants_filt.vcf.gz"
S="/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/baker_633.txt"
REG=[l.split() for l in open(sys.argv[1])]   # name chrom blockstart blockend (bp)
rng=np.random.default_rng(7); rows=[]
for name,c,bs,be in REG:
    bs,be=int(bs),int(be); L=be-bs; fs,fe=max(1,bs-L//3),be+L//3
    cmd=f"{BCF} view -r {c}:{fs}-{fe} -S {S} --force-samples -m2 -M2 -v snps -Ou {V.format(c=c)} 2>/dev/null | {BCF} +fill-tags -Ou -- -t MAF,F_MISSING 2>/dev/null | {BCF} view -i 'INFO/MAF>=0.1 && INFO/F_MISSING<=0.3' -Ou | {BCF} query -f '%POS[\\t%GT]\\n'"
    txt=subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout
    lines=txt.strip().split("\n")
    if len(lines)>1500: lines=[lines[i] for i in sorted(rng.choice(len(lines),1500,replace=False))]
    pos=np.array([int(l.split("\t",1)[0]) for l in lines])
    G=np.array([[{"0/0":0,"0/1":1,"1/0":1,"1/1":2}.get(g,np.nan) for g in l.split("\t")[1:]] for l in lines],float)
    X=G-np.nanmean(G,1,keepdims=True); M=~np.isnan(X); X=np.where(M,X,0)
    cov=X@X.T; n=M.astype(float)@M.T.astype(float); v=np.sqrt(np.maximum(np.diag(cov),1e-12))
    r2=(cov/np.outer(v,v))**2
    np.savez_compressed(f"ld/{name}.npz",pos=pos,r2=r2.astype(np.float32),bs=bs,be=be)
    D=np.abs(pos[:,None]-pos[None,:]); inb=(pos>=bs)&(pos<=be); iu=np.triu_indices(len(pos),1)
    both_in=(inb[:,None]&inb[None,:])[iu]; both_out=(~inb[:,None]&~inb[None,:])[iu]; d=D[iu]; r=r2[iu]
    for lab,sel in (("block",both_in),("flank",both_out)):
        for lo,hi in ((0,1e4),(1e4,1e5),(1e5,1e6),(1e6,1e9)):
            k=sel&(d>=lo)&(d<hi)
            rows.append(dict(region=name,part=lab,dist=f"{int(lo/1e3)}-{'inf' if hi>1e8 else int(hi/1e3)}kb",
                             n_pairs=int(k.sum()),mean_r2=float(r[k].mean()) if k.any() else np.nan,
                             frac_r2_gt_0_2=float((r[k]>0.2).mean()) if k.any() else np.nan))
    print(name,len(pos),"SNPs",flush=True)
pd.DataFrame(rows).to_csv("ld/ld_summary.tsv",sep="\t",index=False); print("LD_DONE")
