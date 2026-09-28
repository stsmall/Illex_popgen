#!/usr/bin/env python
# Clean Ne-through-time plot from Stairway Plot 2 final.summary (baker-633, folded, mu=3e-9, gen=1yr).
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography"
F=f"{D}/stairway_out/Illex_illecebrosus_baker633.final.summary"
yr=[]; med=[]; lo=[]; hi=[]
with open(F) as fh:
    fh.readline()
    for line in fh:
        p=line.rstrip("\n").split("\t")
        if len(p)<9 or p[5]=="" : continue
        try: y=float(p[5]); m=float(p[6]); l=float(p[7]); h=float(p[8])
        except: continue
        yr.append(y); med.append(m); lo.append(l); hi.append(h)
o=np.argsort(yr); yr=np.array(yr)[o]; med=np.array(med)[o]; lo=np.array(lo)[o]; hi=np.array(hi)[o]
# dedup by year
_,idx=np.unique(yr,return_index=True); yr=yr[idx]; med=med[idx]; lo=lo[idx]; hi=hi[idx]
fig,ax=plt.subplots(figsize=(10,6))
ax.fill_between(yr,lo,hi,color="tab:blue",alpha=0.2,label="95% CI")
ax.plot(yr,med,color="tab:blue",lw=1.5,label="Ne median")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Years before present (gen=1yr)"); ax.set_ylabel("Effective population size Ne")
ax.set_title("Illex illecebrosus (baker-633) — Stairway Plot 2\nfolded genome-wide SFS, µ=3e-9, gen=1yr")
ax.legend(); ax.grid(True,which="both",alpha=0.25)
fig.tight_layout(); fig.savefig(f"{D}/stairway_Ne_curve.png",dpi=140,bbox_inches="tight")
print(f"wrote stairway_Ne_curve.png ({len(yr)} time points)")
print(f"Ne range: {med.min():.3e} - {med.max():.3e}")
print(f"most recent (yr={yr[0]:.2f}): Ne={med[0]:.3e}")
print(f"most ancient (yr={yr[-1]:.0f}): Ne={med[-1]:.3e}")
# harmonic-ish summary at decades
for target in [1,10,100,1000,10000,100000,1000000]:
    i=np.argmin(np.abs(yr-target))
    print(f"~{target:>8} yr: Ne_median={med[i]:.3e}  [{lo[i]:.2e}, {hi[i]:.2e}]")
