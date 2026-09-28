import numpy as np
M="/sietch_colab/data_share/illex/alignments/anchorwave_dir/ill_coin_alignment/maf_by_chrom/1.maf"
REG={"BLOCK":(23_000_000,31_000_000),"CONTROL":(33_000_000,41_000_000),"CONTROL2":(70_000_000,78_000_000)}
acc={k:dict(cols=0,aligned=0,mm=0,ref_gap=0,tgt_gap=0) for k in REG}
with open(M) as fh:
    lines=iter(fh)
    for line in lines:
        if not line.startswith("a"): continue
        s1=next(lines).split(); s2=next(lines).split()
        if s1[0]!="s" or s2[0]!="s": continue
        start=int(s1[2]); t1=s1[6].upper(); t2=s2[6].upper()
        # vectorised column walk
        a=np.frombuffer(t1.encode(),dtype="S1"); b=np.frombuffer(t2.encode(),dtype="S1")
        refgap=(a==b"-"); tgtgap=(b==b"-")
        pos=start+np.cumsum(~refgap)-1          # illex coordinate of each column (0-based)
        for k,(lo,hi) in REG.items():
            sel=(pos>=lo)&(pos<hi)&(~refgap)     # columns that carry an illex base in the interval
            al=sel&(~tgtgap)
            acc[k]["cols"]+=int(sel.sum()); acc[k]["aligned"]+=int(al.sum())
            acc[k]["mm"]+=int((al&(a!=b)).sum()); acc[k]["tgt_gap"]+=int((sel&tgtgap).sum())
            # illex-gap columns inside interval (coindetii insertion) -- count where pos in range & refgap
            acc[k]["ref_gap"]+=int(((pos>=lo)&(pos<hi)&refgap).sum())
for k,v in acc.items():
    print(f"{k:9s} illex bp={v['cols']/1e6:.2f}Mb  aligned(both bases)={v['aligned']/1e6:.2f}Mb ({100*v['aligned']/max(v['cols'],1):.1f}%)  "
          f"coindetii-gap={100*v['tgt_gap']/max(v['cols'],1):.1f}%  coin-insertions={v['ref_gap']/1e6:.2f}Mb  "
          f"mismatch/aligned bp = {100*v['mm']/max(v['aligned'],1):.3f}%")
print("DONE_MAF")
