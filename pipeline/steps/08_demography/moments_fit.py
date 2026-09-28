#!/usr/bin/env python
# moments SFS-based demographic inference on the genome-wide folded SFS (baker-633).
# Fit constant / two-epoch / exponential-growth / three-epoch; AIC model comparison;
# convert to Ne + years (mu=3e-9, gen=1yr, L=1.34e9). Repeat with singletons MASKED (drop-singleton test).
import numpy as np, moments
D="/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography"
MU=3e-9; GEN=1.0; L=1340915910          # accessible sites with data (autosomes excl chr2,chr42)
folded=np.loadtxt(f"{D}/gw_folded_sfs.txt")   # bins 0..633 (0=monomorphic), integer counts
N=1266                                    # haploid sample size (633 diploids)

def make_fs(mask_singleton=False):
    data=np.zeros(N+1)
    data[:len(folded)]=folded
    fs=moments.Spectrum(data, mask_corners=False)
    fs.folded=True
    fs.mask[0]=True                       # monomorphic
    fs.mask[N//2+1:]=True                 # folded-away half (634..1266)
    fs.mask[N]=True
    if mask_singleton: fs.mask[1]=True    # drop-singleton test
    return fs

def const(ns): return moments.Demographics1D.snm(ns)
def two_epoch(p,ns): return moments.Demographics1D.two_epoch(p,ns)      # (nu,T)
def growth(p,ns): return moments.Demographics1D.growth(p,ns)           # (nu,T) exp growth 1->nu
def three_epoch(p,ns): return moments.Demographics1D.three_epoch(p,ns) # (nuB,nuF,TB,TF)

MODELS=[("constant",const,None,None,None),
        ("two_epoch",two_epoch,[2,0.1],[1e-3,1e-4],[1e4,10]),
        ("growth",growth,[2,0.1],[1e-3,1e-4],[1e4,10]),
        ("three_epoch",three_epoch,[1,2,0.1,0.05],[1e-3,1e-3,1e-4,1e-4],[1e4,1e4,10,10])]

def fit(fs, label):
    ns=fs.sample_sizes
    print(f"\n===== {label}  (segregating used: {int(fs.sum())}, masked entries: {int(fs.mask.sum())}) =====")
    print(f"{'model':<12} {'k':>2} {'logL':>12} {'AIC':>12}  params->Ne/time")
    results={}
    for name,func,p0,lb,ub in MODELS:
        if p0 is None:
            model=func(ns); k=0
            model=model.fold() if fs.folded else model
            theta=moments.Inference.optimal_sfs_scaling(model,fs)
            ll=moments.Inference.ll_multinom(model,fs)
        else:
            ff=lambda p,n: func(p,n).fold()
            p0p=moments.Misc.perturb_params(p0,fold=0.5,lower_bound=lb,upper_bound=ub)
            popt=moments.Inference.optimize_log(p0p,fs,ff,lower_bound=lb,upper_bound=ub,
                                                maxiter=200,verbose=0)
            model=ff(popt,ns); k=len(popt)
            theta=moments.Inference.optimal_sfs_scaling(model,fs)
            ll=moments.Inference.ll_multinom(model,fs)
        aic=2*k-2*ll
        Nref=theta/(4*MU*L)
        results[name]=dict(ll=ll,aic=aic,theta=theta,Nref=Nref,popt=(None if p0 is None else popt))
        # interpret
        interp=f"Nref={Nref:,.0f}"
        if name in ("two_epoch","growth"):
            nu,T=popt; interp+=f"  nu={nu:.3f}(Ne_cur={nu*Nref:,.0f})  T={T:.4f}(={2*Nref*T*GEN:,.0f} yr)"
        if name=="three_epoch":
            nuB,nuF,TB,TF=popt; interp+=f"  nuB={nuB:.2f} nuF={nuF:.2f}(Ne_cur={nuF*Nref:,.0f}) TB={2*Nref*TB:,.0f}yr TF={2*Nref*TF:,.0f}yr"
        print(f"{name:<12} {k:>2} {ll:>12.1f} {aic:>12.1f}  {interp}")
    best=min(results,key=lambda m:results[m]['aic'])
    print(f"  -> best by AIC: {best} (dAIC vs constant = {results['constant']['aic']-results[best]['aic']:.0f})")
    return results

fs=make_fs(False); r_all=fit(fs,"ALL bins")
fs2=make_fs(True); r_ns=fit(fs2,"SINGLETONS MASKED (drop-singleton test)")
print("\n=== observed vs expected SFS (best growth model, first 12 folded bins, ALL bins fit) ===")
ns=fs.sample_sizes
popt=r_all['growth']['popt']; m=growth(popt,ns).fold(); m=m*r_all['growth']['theta']
print("bin  observed     expected   (ratio)")
for i in range(1,13):
    o=fs[i]; e=m[i]; print(f"{i:>3} {o:>11.0f} {e:>11.0f}   {o/e:.2f}")
