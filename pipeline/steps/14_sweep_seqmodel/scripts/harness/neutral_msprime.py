"""Neutral training sims: msprime under the Illex 'moments' growth demography.

Pure coalescent (no rescaling needed) reproducing genome-wide pi ~ 0.0093 and
Tajima's D ~ -2 -- validated in
``results/slim_validation/msprime_demography_crosscheck.py`` (pi 0.00927, D -1.99),
which is also the reference the SLiM+recapitation sweep flanks were shown to match.

One 2.2 Mb neutral sim yields 11 neutral training windows via
``ms_export.window_offsets`` (the same slicing used for the sweep sims), so the
'neutral' class is produced consistently with the sweep classes' backgrounds.

NB: sweep classes come from Q-rescaled forward SLiM (+recapitate+overlay); this
neutral class is the exact (unrescaled) coalescent. Their pi/Tajima's D match
(validated); any residual rescaling-artifact difference is a Phase-B check.
"""
from __future__ import annotations

import numpy as np
import msprime

# Illex moments growth demography (unscaled; see PLAN_03 Global Constraints)
N0 = 6_808_096       # present size
NREF = 547_928       # ancestral size
TGROW = 769_519      # growth-epoch length (generations)
MU = 3e-9            # per-bp mutation rate
R = 2.1e-9           # per-bp recombination rate


def growth_demography():
    """msprime Demography: present N0 shrinking back at rate g to NREF at TGROW ago."""
    g = np.log(N0 / NREF) / TGROW
    dem = msprime.Demography()
    dem.add_population(name="pop", initial_size=N0, growth_rate=g)
    dem.add_population_parameters_change(time=TGROW, initial_size=NREF, growth_rate=0)
    return dem


def simulate_neutral(seed: int, L: int = 2_200_000, n_dip: int = 350,
                     mu: float = MU, r: float = R):
    """Return a mutated tree sequence: ``2*n_dip`` haplotypes over ``L`` bp.

    Samples are individual-consecutive (haplotypes 2i, 2i+1 = diploid i), matching
    the SLiM reconstruction, so downstream diploid folding (ascertainment) is correct.
    Deterministic from ``seed``.
    """
    dem = growth_demography()
    ts = msprime.sim_ancestry(
        samples={"pop": n_dip}, demography=dem,
        sequence_length=L, recombination_rate=r, random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=mu, random_seed=seed + 1)
    return ts
