import numpy as np
# moments exponential-growth model (present-referenced): N0=present, NREF=ancestral.
NREF=547928.0; N0=6808096.0; T_YEARS=769519.0; GEN=1.0; MU=3e-9
THETA_PER_BP=4*N0*MU          # present-reference scaling
def rho_per_bp(r_per_bp): return 4*N0*r_per_bp
def discoal_demog(N0=N0, NREF=NREF, T_YEARS=T_YEARS, GEN=GEN):
    # discoal NATIVE exponential growth (avoids the -en step-event segfault on small sizes/tiny times):
    #   -g alpha        : present-day exp growth rate; backward in time size = N0*exp(-alpha*tau)
    #   -eg T_scaled 0  : stop growth at the ancestral epoch (size constant = NREF beyond T)
    # time unit = 4*N0 generations; at tau=T_scaled, size = N0*exp(-alpha*T_scaled) = NREF.
    T_scaled=(T_YEARS/GEN)/(4.0*N0)
    alpha=np.log(N0/NREF)/T_scaled
    return ["-g", f"{alpha:.4f}", "-eg", f"{T_scaled:.6f}", "0.0"]
