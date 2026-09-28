#!/usr/bin/env python
# chr2 inversion cline test (Priority 1.3): karyotype frequencies by NAFO division +
# latitude regression. Prediction (panmixia): NO cline -> homogeneous freqs, slope~0.
import pandas as pd, numpy as np
from scipy.stats import chi2_contingency, linregress, pearsonr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

BASE="/sietch_colab/data_share/illex/popgen_data"
D=f"{BASE}/analysis/steps/04_hwe_cline"
k=pd.read_csv(f"{BASE}/analysis/steps/03_karyotype/karyotypes.baker.tsv",sep="\t")
m=pd.read_csv(f"{BASE}/pggpu_illex/popstats/metadata.tsv",sep="\t")[["sample","population","latitude"]]
df=k.merge(m,on="sample").rename(columns={"population":"division"})
df["Bdose"]=df["karyotype"].map({"AA":0,"AB":1,"BB":2})   # B-arrangement allele count
print(f"n={len(df)} baker samples with karyotype+division+latitude", flush=True)

# per-division karyotype counts + arrangement freq + mean lat
tab=df.groupby(["division","karyotype"]).size().unstack(fill_value=0)
for kk in ["AA","AB","BB"]:
    if kk not in tab: tab[kk]=0
tab=tab[["AA","AB","BB"]]
tab["N"]=tab.sum(axis=1)
tab["p_B"]=(2*tab["BB"]+tab["AB"])/(2*tab["N"])
tab["mean_lat"]=df.groupby("division")["latitude"].mean()
tab=tab.sort_values("mean_lat")
print("\n=== karyotype freq by NAFO division (sorted by latitude) ===", flush=True)
print(tab.to_string(), flush=True)

# chi-sq homogeneity of karyotype counts across divisions
chi,p_chi,dof,_=chi2_contingency(tab[["AA","AB","BB"]].values)
print(f"\nChi-sq homogeneity (karyotype x division): chi2={chi:.2f}, df={dof}, p={p_chi:.4g}", flush=True)

# per-individual regression: B-dosage ~ latitude
lr=linregress(df["latitude"],df["Bdose"])
print(f"Per-individual B-dosage ~ latitude: slope={lr.slope:.5f}/deg, r={lr.rvalue:.3f}, p={lr.pvalue:.4g}", flush=True)
# per-division arrangement freq ~ mean latitude
r_div,p_div=pearsonr(tab["mean_lat"],tab["p_B"])
print(f"Per-division p_B ~ mean latitude: r={r_div:.3f}, p={p_div:.4g}", flush=True)

verdict = "NON-CLINAL (panmictic): freqs homogeneous, no latitude effect" if (p_chi>0.05 and lr.pvalue>0.05 and p_div>0.05) else "SOME heterogeneity/cline detected -- inspect"
print(f"\nVERDICT: {verdict}", flush=True)

# plot
fig,ax=plt.subplots(1,2,figsize=(14,5))
bottom=np.zeros(len(tab)); col={"AA":"tab:blue","AB":"tab:orange","BB":"tab:green"}
for kk in ["AA","AB","BB"]:
    fr=tab[kk]/tab["N"]; ax[0].bar(range(len(tab)),fr,bottom=bottom,label=kk,color=col[kk]); bottom+=fr
ax[0].set_xticks(range(len(tab))); ax[0].set_xticklabels([f"{d}\n{la:.0f}N" for d,la in zip(tab.index,tab["mean_lat"])],fontsize=8)
ax[0].set(ylabel="karyotype frequency",title=f"Karyotype freq by NAFO division (chi2 p={p_chi:.2g})"); ax[0].legend()
ax[1].scatter(tab["mean_lat"],tab["p_B"],s=[n/3 for n in tab["N"]])
for d,row in tab.iterrows(): ax[1].annotate(d,(row["mean_lat"],row["p_B"]),fontsize=7)
xs=np.array([tab["mean_lat"].min(),tab["mean_lat"].max()])
ax[1].plot(xs,lr.intercept/1+lr.slope*xs*0+tab["p_B"].mean()+0*xs,":",color="grey")  # flat ref
ax[1].set(xlabel="mean latitude (N)",ylabel="B-arrangement freq",title=f"p_B vs latitude (r={r_div:.2f}, p={p_div:.2g})")
fig.suptitle("chr2 inversion: non-clinal test (baker 633)")
fig.tight_layout(); fig.savefig(f"{D}/cline_test.png",dpi=140,bbox_inches="tight")
tab.to_csv(f"{D}/karyo_by_division.tsv",sep="\t")
print("wrote cline_test.png + karyo_by_division.tsv", flush=True)
