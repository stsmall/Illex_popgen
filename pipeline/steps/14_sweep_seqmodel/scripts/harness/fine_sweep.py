#!/usr/bin/env python3
"""Fine-scale (10kb) pi + Tajima's D across the Z:20.4-21.8Mb sweep zone (illex 330 males),
to locate the valley bottom (= sweep target vicinity) and overlay genes. Also report the
high-frequency-derived-variant count per window (hard-sweep leaves a few near the target)."""
import sys, numpy as np, allel
VCF=sys.argv[1]; A,B,W=20400000,21800000,10000
c=allel.read_vcf(VCF, region=f"Z:{A}-{B}", fields=['variants/POS','calldata/GT'])
pos=c['variants/POS']; gt=allel.GenotypeArray(c['calldata/GT']); ac=gt.count_alleles()
pi,win,_,_=allel.windowed_diversity(pos,ac,size=W,start=A,stop=B)
d,_,_=allel.windowed_tajima_d(pos,ac,size=W,start=A,stop=B)
# high-freq derived (folded proxy: minor-allele freq>=0.7 alt, i.e. alt is common) per window
af=ac.to_frequencies()[:,1]  # ALT freq
hi=(af>=0.5)  # ALT-major sites (post-sweep the beneficial-linked allele is common)
genes=[("unk 20.72-20.74",20717476,20743436),("unk 20.81-20.83",20813328,20829153),
       ("Pasilla 20.83-20.90",20831667,20904350),("Amn 21.47-21.54",21471373,21541796)]
print(f"{'win(Mb)':<10}{'pi':>9}{'TajD':>8}{'nSNP':>6}{'altMaj':>7}  gene")
for i,(a,b) in enumerate(win):
    g=next((nm for nm,gs,ge in genes if gs<=b and ge>=a),"")
    m=(pos>=a)&(pos<=b); nsnp=int(m.sum()); alt=int(hi[m].sum()) if nsnp else 0
    star="  <== pi-min" if np.isfinite(pi[i]) and pi[i]==np.nanmin(pi) else ""
    print(f"{a/1e6:<10.2f}{pi[i]:>9.5f}{d[i]:>8.2f}{nsnp:>6}{alt:>7}  {g}{star}")
mi=int(np.nanargmin(pi))
print(f"\nPI-MIN window: Z:{win[mi][0]/1e6:.2f}-{win[mi][1]/1e6:.2f}Mb  pi={pi[mi]:.5f}  D={d[mi]:.2f}")
di=int(np.nanargmin(d))
print(f"D-MIN  window: Z:{win[di][0]/1e6:.2f}-{win[di][1]/1e6:.2f}Mb  D={d[di]:.2f}  pi={pi[di]:.5f}")
