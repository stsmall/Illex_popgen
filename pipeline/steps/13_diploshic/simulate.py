import sys, numpy as np
sys.path.insert(0,"/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic")
import demography as D
WINLEN=44000; NSUB=11; NHAP=700   # downsampled from baker-633 (discoal n<=~700); 44kb keeps segsites < MAXMUTS
def theta(): return D.THETA_PER_BP*WINLEN
def rho_draw(rng):
    # cM/Mb ~ truncated around ReLERNN 0.21 (0.05..0.6) -> per-bp -> 4N0 r L
    cm=np.clip(rng.normal(0.21,0.08),0.05,0.6); r=cm*1e-8
    return D.rho_per_bp(r)*WINLEN
def base(rng): return [str(NHAP),"1",str(WINLEN),"-t",f"{theta():.3f}","-r",f"{rho_draw(rng):.3f}"]+D.discoal_demog()
def neutral(rng): return base(rng)
def sweep(rng, subwin):
    x=(subwin+0.5)/NSUB                                   # sweep position 0..1
    alpha=2*D.N0*10**rng.uniform(-4,-2)                  # 2Ns via s~1e-4..1e-2
    tau=rng.uniform(0,0.05)                               # fixation time (4N0 units), recent
    args=base(rng)+["-ws",f"{tau:.4f}","-a",f"{alpha:.2f}","-x",f"{x:.4f}"]
    if rng.random()<0.5:                                  # half soft: standing var start freq
        args+=["-f",f"{10**rng.uniform(-3,-1):.5f}"]
    return args, ("soft" if "-f" in args else "hard")
