import sys, os, subprocess
sys.path.insert(0,"/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic")
import demography as D
args=D.discoal_demog()
# native exponential growth: -g <alpha> -eg <T_scaled> 0.0
assert args[0]=="-g" and args[2]=="-eg" and args[4]=="0.0", args
alpha=float(args[1]); Tsc=float(args[3])
# alpha = ln(N0/NREF)/T_scaled ; at tau=T_scaled the size returns to NREF/N0
import math
assert abs(math.exp(-alpha*Tsc) - 547928/6808096) < 1e-4, math.exp(-alpha*Tsc)
assert 0.02 < Tsc < 0.04, Tsc          # 769519/(4*6808096) ~ 0.02826
assert 50 < alpha < 150, alpha         # ln(12.42)/0.02826 ~ 89.2
assert abs(D.THETA_PER_BP - 4*6808096*3e-9) < 1e-9
# SMOKE TEST: discoal (oobfix build) must accept the flags for NEUTRAL, HARD and SOFT sweeps
# at the real downsampled sample size, and produce non-degenerate output.
discoal=os.environ.get("DISCOAL","/home/ssmall/programs/discoal_oobfix/discoal")
def seg(extra):
    r=subprocess.run([discoal,"700","1","50000","-t","4085","-r","2859"]+extra+args,
                     capture_output=True,text=True)
    assert r.returncode==0, f"discoal exit {r.returncode}: {r.stderr[-300:]}"
    s=[int(l.split()[1]) for l in r.stdout.splitlines() if l.startswith("segsites")]
    assert s and s[0]>0, r.stdout[:200]
    return s[0]
n=seg([]); h=seg(["-ws","0.005","-a","500","-x","0.5"]); sf=seg(["-ws","0.005","-a","500","-x","0.5","-f","0.02"])
print(f"demography check PASSED (alpha={alpha:.2f} Tsc={Tsc:.5f}; discoal n700 segsites neutral={n} hard={h} soft={sf})")
