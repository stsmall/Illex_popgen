import os; os.environ.setdefault("MPLCONFIGDIR","/dev/shm/mplcache")
import subprocess, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
BCF="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"
T="/sietch_colab/data_share/illex/popgen_data/analysis/steps/05_karyo_pca_arms/thinned"
S="/sietch_colab/data_share/illex/popgen_data/analysis/steps"
O="/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content"
def gt(c,lo,hi):
    out=subprocess.run([BCF,"query","-r",f"{c}:{int(lo*1e6)}-{int(hi*1e6)}","-f","[%GT\t]\n",f"{T}/chr{c}.thin.vcf.gz"],capture_output=True,text=True).stdout
    return np.array([[np.nan if x[0]=="." else (x[0]!="0")+(x[2]!="0") for x in l.rstrip("\t\n").split("\t")] for l in out.splitlines() if l.strip()],dtype=float)
samples=subprocess.run([BCF,"query","-l",f"{T}/chr1.thin.vcf.gz"],capture_output=True,text=True).stdout.split()
cov=pd.read_csv(f"{S}/12_sex/chrZ_chr1_cov.tsv",sep="\t").set_index("sample")
kar=pd.read_csv(f"{S}/05_karyo_pca_arms/karyo_pca_coords.tsv",sep="\t"); kar=kar[kar.chrom.astype(str)=="2"].set_index("sample")["karyotype"]
# division from per-population bamlists
div={}
import glob
for f in glob.glob(f"{S}/08_demography/per_population/bamlists/*.bamlist.txt"):
    d=os.path.basename(f).split(".")[0]
    for l in open(f): div[os.path.basename(l.strip()).replace(".bam","")]=d
fig,axs=plt.subplots(1,2,figsize=(10,3.6))
for ax,(lab,c,lo,hi) in zip(axs,(("chr1 BLOCK 23-31 Mb","1",23,31),("chr34 BLOCK 29.9-44.7 Mb","34",29.9,44.7))):
    G=gt(c,lo,hi); miss=np.isnan(G).mean(0); het=(G==1).sum(0)/np.maximum((~np.isnan(G)).sum(0),1)
    col=np.nanmean(G,0); col=np.where(np.isnan(col),np.nanmean(G),col); Gi=np.where(np.isnan(G),col,G)
    X=((Gi-col)/np.sqrt(np.clip(col/2*(1-col/2),1e-6,None))).T
    pc1=PCA(n_components=2).fit_transform(X)[:,0]
    g=GaussianMixture(3,n_init=4,random_state=0).fit(pc1[:,None]); lab3=g.predict(pc1[:,None]); order=np.argsort(g.means_.ravel()); rank={o:i for i,o in enumerate(order)}
    cl=np.array([rank[l] for l in lab3])   # 0=majority side ... 2=far side
    df=pd.DataFrame({"sample":samples,"pc1":pc1,"cluster":cl,"het_block":het,"miss_block":miss})
    df["cov_chr1"]=df["sample"].map(cov["chr1_reads"]); df["kar2"]=df["sample"].map(kar); df["div"]=df["sample"].map(div)
    print(f"\n=== {lab} ===  cluster sizes {np.bincount(cl,minlength=3).tolist()}")
    print(df.groupby("cluster").agg(n=("sample","size"),het_block=("het_block","mean"),miss_block=("miss_block","mean"),cov_chr1=("cov_chr1","median")).round(4).to_string())
    print("chr2 karyotype by cluster (fractions):\n"+pd.crosstab(df.cluster,df.kar2,normalize="index").round(2).to_string())
    print("division by cluster (fractions):\n"+pd.crosstab(df.cluster,df["div"],normalize="index").round(2).to_string())
    # HWE-style test: het/hom ratio; for a true inversion middle cluster het >> outer clusters
    h=df.groupby("cluster").het_block.mean(); print(f"het ratio middle/outer: {h.get(1,np.nan)/max(h.get(0,np.nan),h.get(2,np.nan)):.2f}")
    df.to_csv(f"{O}/block_clusters_chr{c}.tsv",sep="\t",index=False)
    ax.hist(pc1,bins=60,color="#0072B2"); ax.set_title(f"{lab}: PC1"); ax.set_xlabel("PC1"); ax.set_ylabel("individuals")
fig.tight_layout(); fig.savefig(f"{O}/block_pc1_hist.png",dpi=150); print("\nwrote block_pc1_hist.png")
