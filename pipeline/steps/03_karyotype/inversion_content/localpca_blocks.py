import subprocess, numpy as np
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
BCF="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"
T="/sietch_colab/data_share/illex/popgen_data/analysis/steps/05_karyo_pca_arms/thinned"
def gt(c, lo, hi):
    out=subprocess.run([BCF,"query","-r",f"{c}:{int(lo*1e6)}-{int(hi*1e6)}","-f","[%GT\t]\n",f"{T}/chr{c}.thin.vcf.gz"],capture_output=True,text=True).stdout
    rows=[[np.nan if x[0]=="." else (x[0]!="0")+(x[2]!="0") for x in line.rstrip("\t\n").split("\t")] for line in out.splitlines() if line.strip()]
    return np.array(rows,dtype=float)
def analyse(lab,c,lo,hi):
    G=gt(c,lo,hi)
    if G.size==0: print(f"{lab}: no SNPs",flush=True); return
    col=np.nanmean(G,0); col=np.where(np.isnan(col),np.nanmean(G),col)
    G=np.where(np.isnan(G),col,G); sd=np.sqrt(np.clip(col/2*(1-col/2),1e-6,None)); X=((G-col)/sd).T
    pcs=PCA(n_components=4).fit(X); Y=pcs.transform(X); ve=pcs.explained_variance_ratio_
    x=Y[:,:1]; bic={k:GaussianMixture(k,n_init=4,random_state=0).fit(x).bic(x) for k in (1,2,3)}
    g3=GaussianMixture(3,n_init=4,random_state=0).fit(x); sep=(g3.means_.max()-g3.means_.min())/np.sqrt(g3.covariances_.ravel().mean())
    cnt=np.bincount(g3.predict(x),minlength=3)[np.argsort(g3.means_.ravel())]; p=(2*cnt[2]+cnt[1])/(2*cnt.sum())
    exp=(np.array([(1-p)**2,2*p*(1-p),p**2])*cnt.sum()).round().astype(int)
    print(f"{lab:26s} chr{c}:{lo}-{hi}Mb SNPs={G.shape[0]:4d} PC1={100*ve[0]:.1f}% PC2={100*ve[1]:.1f}% "
          f"BIC(k1-k3)={bic[1]-bic[3]:+7.1f} (k1-k2={bic[1]-bic[2]:+6.1f}) sep={sep:.1f} 3clust={cnt.tolist()} HWEexp={exp.tolist()}",flush=True)
for lab,c,lo,hi in (("chr2 INVERSION (positive ctl)","2",60.54,79.5),("chr2 collinear","2",10,29),
                    ("chr1 BLOCK","1",23,31),("chr1 control","1",33,41),("chr1 control2","1",70,78),
                    ("chr43 BLOCK","43",45.7,61.6),("chr43 control","43",10,26),
                    ("chr34 BLOCK","34",29.9,44.7),("chr34 control","34",5,20),
                    ("chr38 BLOCK","38",19,31.8),("chr38 control","38",40,53),
                    ("chr21 BLOCK","21",0,9.3),("chr21 control","21",30,40)):
    try: analyse(lab,c,lo,hi)
    except Exception as e: print(lab,"ERR",repr(e)[:120],flush=True)
print("DONE_LOCALPCA")
