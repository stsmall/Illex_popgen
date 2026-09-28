## Illex illecebrosus population-genetic parameters for Bvalcalc
## SMOKE-TEST placeholders except u (illex mutation rate). TAILOR before real use.
## Core parameters
x = 1 # Scaling factor (N,u,r), keep as 1 unless calculating for rescaled simulations
Nanc = 7.75e5 / x # Ancestral Ne (placeholder: pi 0.0093 / (4*u) ~ 7.75e5); tailor
r = 2.086e-9 * x # Fallback crossover rate per bp/gen (genome-avg); overridden per-chunk by --rec_map
u = 3e-9 * x # Mutation rate per bp/gen (illex)
g = 0 * x # Gene conversion initiation rate per bp/gen (set 0 unless a gc_map is supplied)
k = 440 # Gene conversion tract length (bp)
## Basic DFE parameters for ALL sites in annotated regions (Sum must equal 1)
f0 = 0.25 # Proportion effectively neutral 0 <= |2Ns| < 1 (2Ns<5 does not contribute to BGS)
f1 = 0.49 # Proportion weakly deleterious 1 <= |2Ns| < 10
f2 = 0.04 # Proportion moderately deleterious 10 <= |2Ns| < 100
f3 = 0.22 # Proportion strongly deleterious |2Ns| >= 100
## Demography parameters (only used with --pop_change)
Ncur = 2 * Nanc # Current population size
time_of_change = 0.45 * Nanc # Generations ago of Nanc->Ncur change
## Advanced DFE parameters
h = 0.5 # Dominance coefficient
mean, shape, proportion_synonymous = 811/(2*Nanc), 0.347, 0.3 # Gamma DFE [mean s, shape, neutral prop] (--gamma_dfe)
s_breaks = 0, 1/(2*Nanc), 10/(2*Nanc), 100/(2*Nanc), 1 # (--custom_dfe)
bin_proportions = 0.25, 0.25, 0.25, 0.25 # (--custom_dfe)
