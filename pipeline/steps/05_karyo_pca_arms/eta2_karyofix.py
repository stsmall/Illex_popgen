"""Among-karyotype variance (eta^2) of PC1 with the corrected chr2 karyotypes (data/qc/karyotypes_619.tsv, karyotype_final).

eta^2 = between-karyotype sum of squares / total sum of squares of PC1 (one-way ANOVA), over the labelled individuals.
PC1 sign does not matter (eta^2 is invariant to sign and scale), so no orientation is needed.
1) windowed PCA along chr2 (analysis/steps/02_inv_pca/chr2_winpca_coords.npy, 200-kb windows): mean over windows
   inside the breakpoints (60.54-79.50 Mb) and over collinear windows >1 Mb from either breakpoint.
2) per-chromosome PCA of thinned SNPs (analysis/steps/05_karyo_pca_arms/karyo_pca_coords.tsv, karyo_pca_panel.py).
The PCA coordinates are unchanged; only the labels are new. Prints old-label values for comparison.
"""
import subprocess, numpy as np, pandas as pd
A="/sietch_colab/data_share/illex/popgen_data/analysis/steps"
M="/sietch_colab/data_share/illex/popgen_data/illex_manuscript"
old=pd.read_csv(f"{A}/03_karyotype/karyotypes.baker.tsv",sep="\t").set_index("sample")["karyotype"]
new=pd.read_csv(f"{M}/data/qc/karyotypes_619.tsv",sep="\t",comment="#").set_index("sample")["karyotype_final"]
def eta2(x,lab):
    ok=np.isin(lab,["AA","AB","BB"])&np.isfinite(x); x=x[ok]; l=lab[ok]
    sb=sum((l==k).sum()*(x[l==k].mean()-x.mean())**2 for k in ["AA","AB","BB"])
    return sb/((x-x.mean())**2).sum()
# per-chromosome PCA
d=pd.read_csv(f"{A}/05_karyo_pca_arms/karyo_pca_coords.tsv",sep="\t",dtype={"chrom":str})
r=[]
for c,g in d.groupby("chrom"):
    s=g["sample"].to_numpy(); x=g.PC1.to_numpy()
    r.append((c,eta2(x,np.array([old.get(i,"") for i in s])),eta2(x,np.array([new.get(i,"") for i in s]))))
r=pd.DataFrame(r,columns=["chrom","old","new"]).sort_values("new",ascending=False)
print(r.head(4).to_string()); print("max other old/new", r[r.chrom!="2"].old.max(), r[r.chrom!="2"].new.max(), "n chroms", len(r))
# windowed PCA chr2
V=f"{A}/00_callset/filtered/2/variants_filt.vcf.gz"
names=subprocess.run(["/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools","query","-l",V],capture_output=True,text=True).stdout.split()
C=np.load(f"{A}/02_inv_pca/chr2_winpca_coords.npy")[:,:,0]
w=pd.read_csv(f"{A}/02_inv_pca/chr2_winpca_windows.tsv",sep="\t")
lo=np.array([old.get(i,"") for i in names]); ln=np.array([new.get(i,"") for i in names])
eo=np.array([eta2(C[i],lo) for i in range(len(w))]); en=np.array([eta2(C[i],ln) for i in range(len(w))])
inv=(w.start>=60.54e6)&(w.end<=79.50e6); col=(w.end<60.54e6-1e6)|(w.start>79.50e6+1e6)
for nm,e in (("old",eo),("new",en)):
    print(nm,"inv mean/median",np.nanmean(e[inv]).round(4),np.nanmedian(e[inv]).round(4),"col mean/median",np.nanmean(e[col]).round(4),np.nanmedian(e[col]).round(4))
print("windows", int(inv.sum()), int(col.sum()), "n labelled old/new", (lo!="").sum(), (ln!="").sum())
c=w.center
for nm,m in (("overlap",(w.end>60.54e6)&(w.start<79.50e6)),("60-80",(c>=60e6)&(c<=80e6)),("55-85",(c>=55e6)&(c<=85e6)),("center in",(c>=60.54e6)&(c<=79.5e6))):
    print(nm,int(m.sum()),"old",np.nanmean(eo[m]).round(4),"new",np.nanmean(en[m]).round(4))
