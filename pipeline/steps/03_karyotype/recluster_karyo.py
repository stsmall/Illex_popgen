#!/usr/bin/env python
# Corrected chr2 karyotyping: GMM(3) on the saved local-PCA PC1 (inversion_pca.py's
# tertile cut forced a non-HWE 1:1:1 split). Aligns coords rows to VCF sample order,
# writes karyotypes.tsv (sample,PC1,PC2,region_het,karyotype), and a corrected plot.
import numpy as np, pandas as pd, subprocess
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

BASE="/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype"
COORDS=f"{BASE}/chr2_karyo_coords.tsv"
VCF="/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/2/variants_filt.vcf.gz"
BCF="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"

df=pd.read_csv(COORDS, sep="\t")
# VCF sample order (strip mamba-run blank line); align 1:1 to coords rows
samples=[s for s in subprocess.check_output([BCF,"query","-l",VCF]).decode().split("\n") if s.strip()]
assert len(samples)==len(df), f"sample/coords mismatch {len(samples)} vs {len(df)}"
df.insert(0,"sample",samples)

pc1=df["PC1"].values.astype(float).reshape(-1,1)
gm=GaussianMixture(n_components=3, n_init=10, random_state=0).fit(pc1)
lab=gm.predict(pc1); order=np.argsort(gm.means_.ravel())
name={order[0]:"AA", order[1]:"AB", order[2]:"BB"}   # PC1 low->high = AA/AB/BB
df["karyotype"]=[name[l] for l in lab]

het=df["region_het"].values.astype(float)
sizes={k:int((df.karyotype==k).sum()) for k in ["AA","AB","BB"]}
tot=len(df); p=(2*sizes["BB"]+sizes["AB"])/(2*tot)
exp={"AA":(1-p)**2*tot,"AB":2*p*(1-p)*tot,"BB":p*p*tot}
chi=sum((sizes[k]-exp[k])**2/exp[k] for k in sizes)
print("sizes",sizes,"p=%.3f chi2=%.2f"%(p,chi))
for k in ["AA","AB","BB"]:
    print(f"  {k}: n={sizes[k]} het={het[df.karyotype==k].mean():.4f}")

df[["sample","PC1","PC2","region_het","karyotype"]].to_csv(f"{BASE}/karyotypes.tsv", sep="\t", index=False)
print("wrote karyotypes.tsv")

# corrected 3-panel plot
col={"AA":"tab:blue","AB":"tab:orange","BB":"tab:green"}
c=[col[k] for k in df.karyotype]
fig,ax=plt.subplots(1,3,figsize=(20,5.5))
ax[0].scatter(df.PC1,df.PC2,c=c,s=10); ax[0].set(xlabel="PC1",ylabel="PC2",title="chr2:60-80Mb local PCA (GMM karyotype)")
ax[1].hist(df.PC1,bins=60,color="grey"); ax[1].set(xlabel="PC1",ylabel="individuals",title="PC1 (trimodal)")
for k in ["AA","AB","BB"]:
    m=df.karyotype==k; ax[2].scatter(df.PC1[m],het[m],c=col[k],s=10,label=f"{k} n={sizes[k]} H={het[m].mean():.3f}")
ax[2].set(xlabel="PC1",ylabel="region heterozygosity",title="het by karyotype (AB highest)"); ax[2].legend()
fig.suptitle(f"chr2 karyotypes (GMM) — p(B)={p:.3f}, HWE chi2={chi:.2f} (n.s.)  n=653")
fig.tight_layout(); fig.savefig(f"{BASE}/chr2_karyo_corrected.png",dpi=140,bbox_inches="tight")
print("wrote chr2_karyo_corrected.png")
