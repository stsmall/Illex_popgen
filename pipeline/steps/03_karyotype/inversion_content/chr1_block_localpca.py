import subprocess, numpy as np, sys
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
BCF="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"
V="/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/1/variants_filt.vcf.gz"
O="/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content"
def gtmat(region, thin=40, maxsnp=25000):
    p=subprocess.Popen([BCF,"query","-r",region,"-f","[%GT\t]\n",V],stdout=subprocess.PIPE,text=True)
    rows=[]
    for i,line in enumerate(p.stdout):
        if i%thin: continue
        g=line.rstrip("\t\n").split("\t")
        v=np.array([np.nan if x[0]=="." else (x[0]!="0")+(x[2]!="0") for x in g],dtype=float)
        if np.nanmean(v)/2>0.05 and np.nanmean(v)/2<0.95: rows.append(v)   # MAF>5%
        if len(rows)>=maxsnp: break
    p.kill(); return np.array(rows)
out=open(f"{O}/chr1_block_localpca.txt","w")
res={}
for lab,reg in (("BLOCK","1:23000000-31000000"),("CONTROL","1:33000000-41000000")):
    G=gtmat(reg); n_snp,n_ind=G.shape
    col=np.nanmean(G,0); G=np.where(np.isnan(G),col,G); G=(G-col)/np.sqrt(col/2*(1-col/2)+1e-9)
    pcs=PCA(n_components=5).fit(G.T); X=pcs.transform(G.T)
    ve=pcs.explained_variance_ratio_
    bic={k:GaussianMixture(k,n_init=5,random_state=0).fit(X[:,:1]).bic(X[:,:1]) for k in (1,2,3,4)}
    best=min(bic,key=bic.get)
    gm=GaussianMixture(3,n_init=5,random_state=0).fit(X[:,:1]); lab3=gm.predict(X[:,:1])
    cnt=np.bincount(lab3,minlength=3); means=gm.means_.ravel(); order=np.argsort(means)
    # for a 3-cluster inversion: middle cluster ~ heterozygotes; check HWE-like proportions
    c=cnt[order]; p=(2*c[2]+c[1])/(2*c.sum())
    exp=np.array([(1-p)**2,2*p*(1-p),p**2])*c.sum()
    print(f"{lab} {reg}: SNPs={n_snp} ind={n_ind} | PC1 var={100*ve[0]:.1f}% PC2={100*ve[1]:.1f}% PC3={100*ve[2]:.1f}% "
          f"| GMM-BIC best k={best} (BIC k1={bic[1]:.0f} k2={bic[2]:.0f} k3={bic[3]:.0f} k4={bic[4]:.0f}) "
          f"| 3-cluster sizes (sorted by PC1) {c.tolist()} vs HWE expectation {exp.round(0).astype(int).tolist()} (p_alt={p:.3f})", file=out, flush=True)
    np.save(f"{O}/chr1_{lab}_pc.npy", X)
print("DONE", file=out); out.close()
