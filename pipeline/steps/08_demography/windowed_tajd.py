#!/usr/bin/env python
# 3.1 Genome-wide windowed diversity (50kb windows, 10kb step) from per-chrom thetaStat win50k.pestPG.
# Manhattan-style plot of Tajima's D + pi across the genome (baker-633, folded, accessible sites).
import glob, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/thetas"
CHRS=[str(c) for c in range(1,46)]+["Z"]
cols=["paren","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites"]
frames=[]
for c in CHRS:
    f=f"{D}/baker.{c}.win50k.pestPG"
    if not os.path.exists(f): continue
    df=pd.read_csv(f,sep="\t",header=0,names=cols)
    df=df[df["nSites"]>=5000].copy()          # require >=5kb data in the 50kb window
    df["chrom"]=c; df["pi"]=df["tP"]/df["nSites"]; df["thetaW"]=df["tW"]/df["nSites"]
    frames.append(df[["chrom","WinCenter","pi","thetaW","Tajima","nSites"]])
g=pd.concat(frames,ignore_index=True)
g.to_csv(f"{D.replace('/thetas','')}/windowed_diversity.tsv",sep="\t",index=False)
print(f"{len(g)} windows across {g['chrom'].nunique()} chroms (>=5kb data)")
print(f"Tajima's D: mean={g['Tajima'].mean():.3f} median={g['Tajima'].median():.3f} sd={g['Tajima'].std():.3f}")
print(f"pi: mean={g['pi'].mean():.6f}  frac windows D<0: {(g['Tajima']<0).mean():.3f}")

# cumulative genome coordinate for Manhattan plot
order=CHRS
g["chrom"]=pd.Categorical(g["chrom"],categories=order,ordered=True)
g=g.sort_values(["chrom","WinCenter"])
offset=0; xs=[]; ticks=[]; ticklab=[]
for c in order:
    sub=g[g["chrom"]==c]
    if len(sub)==0: continue
    xs.append(sub["WinCenter"].values+offset)
    ticks.append(offset+sub["WinCenter"].max()/2); ticklab.append(c)
    offset+=sub["WinCenter"].max()+1e6
g["gx"]=np.concatenate(xs)

fig,ax=plt.subplots(2,1,figsize=(18,8),sharex=True)
colmap={c:("#3b6ea5" if i%2==0 else "#a5c8e1") for i,c in enumerate(order)}
cols_pts=g["chrom"].map(colmap)
ax[0].scatter(g["gx"],g["Tajima"],s=1.5,c=cols_pts,rasterized=True)
ax[0].axhline(0,color="grey",lw=0.6); ax[0].axhline(g["Tajima"].mean(),color="red",lw=0.8,ls="--",label=f"mean={g['Tajima'].mean():.2f}")
ax[0].set(ylabel="Tajima's D"); ax[0].legend(loc="lower right",fontsize=9)
ax[0].set_title("Illex illecebrosus (baker-633) genome-wide windowed diversity — 50kb win / 10kb step, folded, accessible")
ax[1].scatter(g["gx"],g["pi"],s=1.5,c=cols_pts,rasterized=True)
ax[1].axhline(g["pi"].mean(),color="red",lw=0.8,ls="--",label=f"mean pi={g['pi'].mean():.4f}")
ax[1].set(ylabel="pi (per site)"); ax[1].legend(loc="upper right",fontsize=9)
ax[1].set_xticks(ticks); ax[1].set_xticklabels(ticklab,fontsize=7); ax[1].set_xlabel("chromosome")
fig.tight_layout(); out=f"{D.replace('/thetas','')}/windowed_diversity.png"
fig.savefig(out,dpi=130,bbox_inches="tight"); print(f"wrote {out}")
