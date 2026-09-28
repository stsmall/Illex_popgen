#!/usr/bin/env python
# Per-chromosome PCA colored by chr2 karyotype (chr2-specificity test).
# If the inversion drives the structure, only chr2 shows 3 karyotype clusters;
# other chromosomes show samples mixed. Quantify via variance-explained of PC1 by karyotype.
import os, glob, subprocess, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pg_gpu import HaplotypeMatrix, decomposition
from scipy.stats import f_oneway

D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/05_karyo_pca_arms"
BCF="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"
karyo=pd.read_csv("/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/karyotypes.baker.tsv",sep="\t")
kmap=dict(zip(karyo["sample"],karyo["karyotype"]))
col={"AA":"tab:blue","AB":"tab:orange","BB":"tab:green"}

def chrom_key(f):
    c=os.path.basename(f).replace("chr","").replace(".thin.vcf.gz","")
    return (0,int(c)) if c.isdigit() else (1,c)   # numeric first, then Z
vcfs=sorted(glob.glob(f"{D}/thinned/chr*.thin.vcf.gz"), key=chrom_key)
print(f"{len(vcfs)} thinned chrom VCFs", flush=True)

rows=[]; results={}
for vcf in vcfs:
    c=os.path.basename(vcf).replace("chr","").replace(".thin.vcf.gz","")
    try:
        h=HaplotypeMatrix.from_vcf(vcf)
        samples=[s for s in subprocess.check_output([BCF,"query","-l",vcf]).decode().split("\n") if s.strip()]
        nind=len(samples)
        coords,ev=decomposition.randomized_pca(h,n_components=4)
        coords=coords.get() if hasattr(coords,"get") else coords
        if coords.shape[0]==2*nind: coords=(coords[:nind]+coords[nind:2*nind])/2.0  # fold hap pairs -> per-individual
        kt=np.array([kmap.get(s,"NA") for s in samples])
        pc1,pc2=coords[:,0],coords[:,1]
        # variance of PC1 explained by karyotype (ANOVA eta^2)
        grp=[pc1[kt==k] for k in ["AA","AB","BB"] if (kt==k).sum()>1]
        if len(grp)==3:
            F,p=f_oneway(*grp)
            ss_between=sum(len(x)*(x.mean()-pc1.mean())**2 for x in grp); ss_tot=((pc1-pc1.mean())**2).sum()
            eta2=ss_between/ss_tot if ss_tot>0 else 0
        else: F,p,eta2=(np.nan,np.nan,np.nan)
        results[c]=dict(n=h.num_variants,eta2=eta2,F=F,p=p)
        for s,a,b,k in zip(samples,pc1,pc2,kt): rows.append((c,s,a,b,k))
        print(f"chr{c}: {h.num_variants} SNPs, PC1~karyotype eta^2={eta2:.3f} (p={p:.1e})",flush=True)
    except Exception as e:
        print(f"chr{c}: FAILED {e}",flush=True)

df=pd.DataFrame(rows,columns=["chrom","sample","PC1","PC2","karyotype"])
df.to_csv(f"{D}/karyo_pca_coords.tsv",sep="\t",index=False)

# panel
chroms=[c for c in df.chrom.unique()]
n=len(chroms); ncol=7; nrow=int(np.ceil(n/ncol))
fig,ax=plt.subplots(nrow,ncol,figsize=(3*ncol,3*nrow))
ax=np.array(ax).ravel()
for i,c in enumerate(chroms):
    d=df[df.chrom==c]
    ax[i].scatter(d.PC1,d.PC2,c=[col.get(k,"grey") for k in d.karyotype],s=4)
    e=results.get(c,{}).get("eta2",np.nan)
    ax[i].set_title(f"chr{c}  eta2={e:.2f}",fontsize=8,color=("red" if (e==e and e>0.3) else "black"))
    ax[i].set_xticks([]); ax[i].set_yticks([])
for j in range(n,len(ax)): ax[j].axis("off")
fig.suptitle("Per-chromosome PCA colored by chr2 karyotype (AA=blue AB=orange BB=green) — only chr2 clusters",fontsize=12)
fig.tight_layout(); fig.savefig(f"{D}/karyo_pca_panel.png",dpi=120,bbox_inches="tight")
print("wrote karyo_pca_panel.png + karyo_pca_coords.tsv",flush=True)
# summary table
s=pd.DataFrame(results).T.sort_values("eta2",ascending=False)
print("\nPC1-variance-explained-by-karyotype (eta^2), top:\n", s[["n","eta2","p"]].head(8).to_string())
