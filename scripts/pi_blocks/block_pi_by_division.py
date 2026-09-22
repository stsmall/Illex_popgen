"""Are the 22 elevated-pi blocks present in every sampling division, or driven by a subset?
Block calls reproduce the genome-wide screen (pooled 10-kb pi, rolling median > 2x chromosome median, >= 1 Mb).
For each block and division: pi(block)/pi(flanks, 5 Mb each side) from the per-division ANGSD 50-kb windows."""
import pandas as pd, numpy as np
W="/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/bgs_diagnostic_windows.tsv"
P="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/per_population/tables/perpop_windows_expanded.tsv"
w=pd.read_csv(W,sep="\t"); w["chrom"]=w.chrom.astype(str); w["mid"]=(w.start+w.stop)/2
blocks=[]
for c,g in w.groupby("chrom"):
    g=g.sort_values("mid"); med=g.pi.median(); r=g.pi.rolling(61,center=True,min_periods=20).median()
    hi=(r>2*med).to_numpy(); mids=g.mid.to_numpy(); i=0
    while i<len(hi):
        if hi[i]:
            j=i
            while j+1<len(hi) and hi[j+1]: j+=1
            if mids[j]-mids[i]>=1e6: blocks.append((c,mids[i],mids[j]))
            i=j+1
        else: i+=1
b=pd.DataFrame(blocks,columns=["chrom","start","end"]); print(len(b),"blocks")
p=pd.read_csv(P,sep="\t"); p["chrom"]=p.chrom.astype(str)
rows=[]
for bl in b.itertuples():
    q=p[p.chrom==bl.chrom]
    ins=q[(q.WinCenter>=bl.start)&(q.WinCenter<=bl.end)]
    fl=q[((q.WinCenter>=bl.start-5e6)&(q.WinCenter<bl.start))|((q.WinCenter>bl.end)&(q.WinCenter<=bl.end+5e6))]
    for d in sorted(q.division.unique()):
        a=ins[ins.division==d]; f=fl[fl.division==d]
        rows.append(dict(block=f"chr{bl.chrom}:{bl.start/1e6:.1f}-{bl.end/1e6:.1f}",division=d,
                         pi_block=a.pi.mean(),pi_flank=f.pi.mean(),ratio=a.pi.mean()/f.pi.mean(),
                         D_block=a.tajd.mean(),D_flank=f.tajd.mean()))
r=pd.DataFrame(rows); r.to_csv("block_pi_by_division.tsv",sep="\t",index=False)
s=r.groupby("block").agg(ratio_min=("ratio","min"),ratio_median=("ratio","median"),ratio_max=("ratio","max"),
    cv=("ratio",lambda x:x.std()/x.mean()),n_div_gt1_5=("ratio",lambda x:(x>1.5).sum()),n_div=("ratio","size"))
print(s.round(2).to_string())
c1=r[r.block.str.startswith("chr1:")]; print("\nchr1 block by division:\n", c1[["division","pi_block","pi_flank","ratio","D_block","D_flank"]].round(4).to_string(index=False))
