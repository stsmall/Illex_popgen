#!/usr/bin/env python
# 3.3 Pool per-chrom folded SFS -> genome-wide folded SFS for Stairway Plot 2.
# baker-633 (1266 haploid), folded, accessible sites. POOL = autosomes excl chr2 (inversion) + chr42 (sex).
# realSFS -fold output = 1267 floats; folded mass in bins 0..633 (0=invariant), 634..1266 = 0.
import numpy as np, os
D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/thetas"
N_DIP=633; N_HAP=2*N_DIP                      # 1266 chromosomes
NFOLD=N_DIP+1                                 # 634 folded categories (0..633)
POOL=[str(c) for c in range(1,46) if c not in (2,42)]   # autosomes excl inversion+sex
tot=np.zeros(1267)
used=[]
for c in POOL:
    f=f"{D}/baker.{c}.sfs"
    if not os.path.exists(f): continue
    v=np.array(open(f).read().split(),dtype=float)
    if v.size!=1267: print(f"WARN chr{c}: {v.size} entries"); continue
    tot+=v; used.append(c)
folded=tot[:NFOLD].copy()                     # bins 0..633
L=int(round(tot.sum()))                        # total sites (invariant + variant) with data
seg=int(round(folded[1:].sum()))               # segregating sites (exclude bin0)
print(f"pooled {len(used)} chroms (autosomes excl chr2,chr42)")
print(f"total sites L = {L:,}")
print(f"invariant (bin0) = {int(round(folded[0])):,}  ({folded[0]/L*100:.2f}%)")
print(f"segregating = {seg:,}")
print(f"folded SFS bins 1..10: {[int(round(x)) for x in folded[1:11]]}")
print(f"singletons/doubletons ratio = {folded[1]/folded[2]:.2f} (>1.5 = expansion-like excess of rare)")
# save pooled folded SFS (integer counts, bins 0..633)
np.savetxt(f"{os.path.dirname(D)}/gw_folded_sfs.txt",
           np.round(folded).astype(int)[None,:],fmt="%d",delimiter="\t")
# Stairway-ready: bins 1..633 (drop monomorphic), integer
with open(f"{os.path.dirname(D)}/gw_folded_sfs.stairway.txt","w") as fh:
    fh.write("\t".join(str(int(round(x))) for x in folded[1:NFOLD])+"\n")
print(f"\nwrote gw_folded_sfs.txt (bins 0..633) + gw_folded_sfs.stairway.txt (bins 1..633)")
print(f"For Stairway blueprint: nseq={N_HAP}, L={L}, {NFOLD-1} folded bins")
