#!/usr/bin/env python
# Depth-control recluster: local PCA + GMM karyotyping on the DEPTH-EQUALIZED (~1x)
# downsampled chr2:60-80Mb genotypes. Compare to full-depth karyotypes.baker.tsv.
# If clusters persist + calls concordant -> inversion is NOT a depth artifact.
import numpy as np, pandas as pd, subprocess
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from pg_gpu import HaplotypeMatrix, decomposition

BASE="/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype"
DSVCF=f"{BASE}/downsample/chr2inv_ds.vcf.gz"
FULL=f"{BASE}/karyotypes.baker.tsv"
BCF="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"
REGION="2:60000000-80000000"

h=HaplotypeMatrix.from_vcf(DSVCF, region=REGION)
samples=[s for s in subprocess.check_output([BCF,"query","-l",DSVCF]).decode().split("\n") if s.strip()]
# normalize sample names (strip path/.ds.bam if mpileup used filenames)
samples=[s.split("/")[-1].replace(".ds.bam","").replace(".bam","") for s in samples]
n=h.num_haplotypes//2
print(f"[ds] {h.num_variants} SNPs, {n} individuals", flush=True)

coords,ev=decomposition.randomized_pca(h, n_components=6)
coords=coords.get() if hasattr(coords,"get") else coords
# fold haplotype pairs -> per-individual (pg_gpu [hap0_all, hap1_all] layout)
if coords.shape[0]==2*n:
    coords=(coords[:n]+coords[n:2*n])/2.0
pc1=coords[:,0]; pc2=coords[:,1]

# per-individual region heterozygosity from downsampled GTs
gt=subprocess.check_output([BCF,"query","-f","[%GT ]\n",DSVCF]).decode().strip().split("\n")
het=np.zeros(n)
tot=np.zeros(n)
for line in gt:
    a=line.split()
    for i,g in enumerate(a):
        if g in ("0/0","1/1","0|0","1|1"): tot[i]+=1
        elif g in ("0/1","1/0","0|1","1|0"): tot[i]+=1; het[i]+=1
reg_het=np.where(tot>0, het/tot, 0.0)

gm=GaussianMixture(n_components=3,n_init=10,random_state=0).fit(pc1.reshape(-1,1))
lab=gm.predict(pc1.reshape(-1,1)); order=np.argsort(gm.means_.ravel())
name={order[0]:"AA",order[1]:"AB",order[2]:"BB"}
calls=np.array([name[l] for l in lab])
ds=pd.DataFrame({"sample":samples,"PC1":pc1,"PC2":pc2,"region_het":reg_het,"karyotype_ds":calls})
ds.to_csv(f"{BASE}/karyotypes_downsampled.tsv",sep="\t",index=False)

import collections; c=collections.Counter(calls); t=len(calls)
p=(2*c["BB"]+c["AB"])/(2*t)
exp={"AA":(1-p)**2*t,"AB":2*p*(1-p)*t,"BB":p*p*t}
chi=sum((c[k]-exp[k])**2/exp[k] for k in c)
print(f"[ds] sizes AA={c['AA']} AB={c['AB']} BB={c['BB']}  p(B)={p:.3f}  HWE chi2={chi:.2f}", flush=True)
for k in ["AA","AB","BB"]:
    m=calls==k; print(f"   {k}: n={c[k]} het={reg_het[m].mean():.4f}", flush=True)

# concordance vs full-depth
full=pd.read_csv(FULL,sep="\t")[["sample","karyotype"]]
m=ds.merge(full,on="sample")
# karyotype labels may be flipped (AA<->BB) since PC1 sign is arbitrary; align by best match
def concord(a,b): return (a==b).mean()
flip={"AA":"BB","BB":"AA","AB":"AB"}
c1=concord(m.karyotype_ds,m.karyotype)
c2=concord(m.karyotype_ds.map(flip),m.karyotype)
best=max(c1,c2); m["ds_aligned"]=m.karyotype_ds if c1>=c2 else m.karyotype_ds.map(flip)
print(f"[ds] concordance with full-depth calls = {best*100:.1f}% (n={len(m)})", flush=True)
ct=pd.crosstab(m["ds_aligned"],m["karyotype"]); print(ct, flush=True)

col={"AA":"tab:blue","AB":"tab:orange","BB":"tab:green"}
fig,ax=plt.subplots(1,3,figsize=(20,5.5))
ax[0].scatter(pc1,pc2,c=[col[k] for k in calls],s=12); ax[0].set(title="downsampled ~1x local PCA (GMM karyotype)",xlabel="PC1",ylabel="PC2")
ax[1].hist(pc1,bins=60,color="grey"); ax[1].set(title=f"PC1 (uniform depth) p(B)={p:.3f} HWE chi2={chi:.2f}",xlabel="PC1")
for k in ["AA","AB","BB"]:
    mm=calls==k; ax[2].scatter(pc1[mm],reg_het[mm],c=col[k],s=12,label=f"{k} n={c[k]} H={reg_het[mm].mean():.3f}")
ax[2].set(title=f"het by karyotype (concordance vs full-depth {best*100:.0f}%)",xlabel="PC1",ylabel="region het"); ax[2].legend()
fig.suptitle("chr2 karyotyping at UNIFORM ~1x depth (depth control) — baker 633")
fig.tight_layout(); fig.savefig(f"{BASE}/chr2_karyo_downsampled.png",dpi=140,bbox_inches="tight")
print("wrote karyotypes_downsampled.tsv + chr2_karyo_downsampled.png", flush=True)
