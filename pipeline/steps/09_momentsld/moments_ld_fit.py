#!/usr/bin/env python
# momentsLD demographic fit from pg_gpu LD statistics (33 x 2Mb regions, baker-633, flat 1cM/Mb).
# Fit constant(snm) / two-epoch / exponential-growth; AIC compare. Ne is LAST param (rho=4Ne*r).
# NOTE: absolute Ne scales with the assumed flat recombination rate (1e-8/bp); shape/expansion robust.
import pickle, numpy as np, moments.LD
MD="/sietch_colab/data_share/illex/popgen_data/analysis/steps/09_momentsld"
d=pickle.load(open(f"{MD}/ld_regions.pkl","rb"))
regions=d["regions"]; r_bins=np.array(d["r_bins"]); MU=d["mu"]; GEN=1.0
print(f"regions={len(regions)}  r_bins={list(r_bins)}")
mv=moments.LD.Parsing.bootstrap_data(regions)
data=[mv["means"], mv["varcovs"]]
print(f"means: {len(mv['means'])} bin-vectors; stats={mv.get('stats')}")

def run(name, func, p0, lb, ub):
    try:
        opt,LL=moments.LD.Inference.optimize_log_fmin(
            p0, data, [func], rs=r_bins, lower_bound=lb, upper_bound=ub,
            verbose=0, maxiter=300)
        k=len(p0); aic=2*k-2*LL
        return dict(opt=opt,LL=LL,k=k,aic=aic)
    except Exception as e:
        print(f"  {name} FAILED: {type(e).__name__}: {e}")
        return None

MODELS=[
 ("constant", moments.LD.Demographics1D.snm,      [1e4],           [1e2],        [1e9]),
 ("two_epoch",moments.LD.Demographics1D.two_epoch,[2,0.1,1e4],     [1e-3,1e-4,1e2],[1e4,10,1e9]),
 ("growth",   moments.LD.Demographics1D.growth,   [2,0.1,1e4],     [1e-3,1e-4,1e2],[1e4,10,1e9]),
]
res={}
print(f"\n{'model':<11} {'k':>2} {'logL':>12} {'AIC':>12}  interpretation")
for name,func,p0,lb,ub in MODELS:
    r=run(name,func,p0,lb,ub)
    if r is None: continue
    res[name]=r; o=r["opt"]; Ne=o[-1]
    interp=f"Ne_anc={Ne:,.0f}"
    if name in ("two_epoch","growth"):
        nu,T=o[0],o[1]; interp+=f"  nu={nu:.3f} (Ne_cur={nu*Ne:,.0f})  T={T:.4f} (={2*Ne*T*GEN:,.0f} yr)"
    print(f"{name:<11} {r['k']:>2} {r['LL']:>12.1f} {r['aic']:>12.1f}  {interp}")
if res:
    best=min(res,key=lambda m:res[m]['aic'])
    dref=res.get('constant',{}).get('aic',float('nan'))
    print(f"\n-> best by AIC: {best}  (dAIC vs constant = {dref-res[best]['aic']:.1f})")
    import json
    json.dump({k:{'opt':list(map(float,v['opt'])),'LL':float(v['LL']),'aic':float(v['aic'])} for k,v in res.items()},
              open(f"{MD}/momentsld_fit.json","w"), indent=2)
    print("wrote momentsld_fit.json")
