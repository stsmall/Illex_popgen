#!/usr/bin/env python
# Assign sex from chrZ coverage ratio (GMM 2-cluster). High Z_ratio=homogametic(ZZ,Z-diploid);
# low=heterogametic(ZW/Z0,Z-haploid). Under Z nomenclature: ZZ=male, heterogametic=female.
import numpy as np, pandas as pd
from sklearn.mixture import GaussianMixture
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SD="/sietch_colab/data_share/illex/popgen_data/analysis/steps/12_sex"
df=pd.read_csv(f"{SD}/sex_cov.tsv",sep="\t")
x=df["Z_ratio"].values.reshape(-1,1)
gm=GaussianMixture(n_components=2,n_init=10,random_state=0).fit(x)
lab=gm.predict(x); means=gm.means_.flatten()
hi=int(np.argmax(means))                      # high-Z cluster = homogametic ZZ
df["cluster"]=lab
df["karyotype_Z"]=np.where(lab==hi,"ZZ","ZW")  # ZZ=homogametic Z-diploid; ZW=heterogametic Z-haploid
df["sex"]=np.where(lab==hi,"male","female")    # ZZ=male, ZW=female (Z-system inference)
# guard: intermediate/ambiguous (posterior<0.9) flagged
post=gm.predict_proba(x).max(axis=1); df["confident"]=post>=0.9
df.to_csv(f"{SD}/sex_assignment.tsv",sep="\t",index=False)
nZZ=(df.karyotype_Z=="ZZ").sum(); nZW=(df.karyotype_Z=="ZW").sum()
print(f"cluster means (Z_ratio): {sorted(means.round(3))}")
print(f"ZZ (homogametic, male, Z-diploid): {nZZ}")
print(f"ZW (heterogametic, female, Z-haploid): {nZW}")
print(f"ambiguous (posterior<0.9): {(~df.confident).sum()}")
print(f"mean Z_ratio ZZ={df[df.karyotype_Z=='ZZ'].Z_ratio.mean():.3f}  ZW={df[df.karyotype_Z=='ZW'].Z_ratio.mean():.3f}  (ratio {df[df.karyotype_Z=='ZZ'].Z_ratio.mean()/df[df.karyotype_Z=='ZW'].Z_ratio.mean():.2f})")
fig,ax=plt.subplots(1,2,figsize=(12,4.5))
for k,c in [("ZZ","tab:blue"),("ZW","tab:red")]:
    ax[0].hist(df[df.karyotype_Z==k].Z_ratio,bins=40,alpha=0.6,label=f"{k} (n={ (df.karyotype_Z==k).sum() })",color=c)
ax[0].axvline((means.min()+means.max())/2,ls="--",color="grey"); ax[0].set(xlabel="chrZ / autosome coverage ratio",ylabel="individuals",title="Sex from chrZ coverage (GMM)"); ax[0].legend()
ax[1].scatter(df.Z_ratio,df.chr42_ratio,c=(df.karyotype_Z=="ZZ").map({True:"tab:blue",False:"tab:red"}),s=8)
ax[1].set(xlabel="chrZ ratio",ylabel="chr42 ratio",title="chrZ (bimodal=sex) vs chr42 (unimodal=autosomal)")
fig.tight_layout(); fig.savefig(f"{SD}/sex_assignment.png",dpi=140)
print("wrote sex_assignment.tsv + sex_assignment.png")
