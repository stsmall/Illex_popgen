#!/usr/bin/env python
"""Coalescent sanity check of the Illex moments growth demography (no sweep, no rescaling).
Confirms the demographic parameters reproduce pi ~ 0.0093 and Tajima's D ~ -2.
This gives the reference targets the SLiM+recapitation far-windows must match."""
import msprime
import numpy as np

N0 = 6_808_096      # present
NREF = 547_928      # ancestral
TGROW = 769_519     # growth epoch length (gens)
MU = 3e-9
R = 2.1e-9
L = 1_100_000
n_dip = 350         # 700 haplotypes
SEED = 42

g = np.log(N0 / NREF) / TGROW
print(f"growth rate g = {g:.6e} per gen")

dem = msprime.Demography()
# Present population size N0, exponentially shrinking back in time at rate g,
# until TGROW ago it reaches NREF, then constant NREF.
dem.add_population(name="pop", initial_size=N0, growth_rate=g)
dem.add_population_parameters_change(time=TGROW, initial_size=NREF, growth_rate=0)

ts = msprime.sim_ancestry(
    samples={"pop": n_dip},
    demography=dem,
    sequence_length=L,
    recombination_rate=R,
    random_seed=SEED,
)
ts = msprime.sim_mutations(ts, rate=MU, random_seed=SEED)
print("n samples (haplotypes):", ts.num_samples)
print("n sites:", ts.num_sites)

# genome-wide pi and Tajima's D
pi = ts.diversity(mode="site")
taj = ts.Tajimas_D(mode="site")
print(f"genome-wide pi   = {pi:.6f}   (target ~0.0093)")
print(f"genome-wide TajD = {taj:.4f}   (target ~-2.0)")

# windowed (10 windows) to see baseline flatness
wins = np.linspace(0, L, 11)
piw = ts.diversity(windows=wins, mode="site")
tajw = ts.Tajimas_D(windows=wins, mode="site")
print("windowed pi:", np.round(piw, 5))
print("windowed TajD:", np.round(tajw, 3))
