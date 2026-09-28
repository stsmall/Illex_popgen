#!/usr/bin/env python3
"""Design-informed HMM/HSMM params (no fresh calibration sims).

The fresh-sim calibration is infeasible on this box (dense-ARG genotype_matrix
>2.7 h/segment). Per user decision (2026-09-15), build params.json from validated
constants instead:

  emission  : prior_corrected posteriors -- emission_logprob divides the collapsed
              diploSHIC posterior [N,L,C] by the balanced 5-class training prior
              (0.2,0.4,0.4), i.e. Bayes-inverts posterior->likelihood. The illexModel
              confusion structure is already carried in the posteriors; no separate
              matrix needed. This is the module default and is exactly what a fit would
              have chosen (prior_corrected wins emission_accuracy on any balanced set).
  dwell     : d_C / d_L from the ~400 kb diploSHIC sweep footprint (winSize 1.1 Mb /
              11x100 kb; memory diploshic-scan-design) x the genome length-weighted mean
              ReLERNN rate MEAN_RATE=2.51859e-9 M/bp. Center core ~150 kb half-width;
              linked tail ~550 kb (half the classification window).
  recomb    : per-window rate from the ReLERNN sex-avg map (config/mapdf.tsv) inside
              decode; MEAN_RATE fallback for chr2/chr42 (absent from the map).
  null      : filled later by build_empirical_null.py (within-chrom block permutation
              of the empirical scan -> genome-max distribution -> FDR threshold).

Reuses fit_params.build_params() so the JSON schema is identical to a real fit.
"""
import sys, json, os
D = "/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel"
sys.path.insert(0, f"{D}/calibration")
sys.path.insert(0, f"{D}/scripts")
import fit_params as fp   # noqa: E402

MEAN_RATE = fp.MEAN_RATE                       # 2.51859e-9 M/bp

# ---- design dwell scales, from the ~400 kb footprint ----
CENTER_HALF_BP = 150_000                        # hard/soft core half-width
LINKED_HALF_BP = 550_000                        # linked tail half-width (~half the 1.1 Mb window)
d_C = CENTER_HALF_BP * MEAN_RATE               # Morgans
d_L = LINKED_HALF_BP * MEAN_RATE
d_L = max(d_L, d_C * 1.5)

hsmm_dur = {"LL": (d_L - d_C, 0.6), "C": (d_C, 0.6), "LR": (d_L - d_C, 0.6)}
# weakly-informative topology priors (no sims): modest direct-entry / abort rates,
# clipped by build_params to its safe ranges.
topo = {"narrow_entry": 0.10, "abort_left": 0.15, "abort_center": 0.15}
MAX_DUR_WINDOWS = 11                             # footprint ~1.1 Mb / 100 kb windows

params = fp.build_params(d_C, d_L, "prior_corrected", hsmm_dur, topo,
                         MAX_DUR_WINDOWS, mean_win_gap=100_000)
params["null"] = {}                             # thresholds added by build_empirical_null.py
params["labels"]["provenance"] = (
    "DESIGN-INFORMED (no fresh calibration sims; fresh-sim scan infeasible on this box). "
    f"d_C/d_L from ~400 kb diploSHIC footprint x mean ReLERNN rate {MEAN_RATE:.5g} M/bp; "
    "emission=prior_corrected posteriors (balanced 5-class prior 0.2/0.4/0.4); "
    "recomb=ReLERNN sex-avg map (chr2/chr42 constant-rate fallback); "
    "null=empirical within-chrom block permutation FDR.")

out = f"{D}/calibration/params.json"
json.dump(params, open(out, "w"), indent=2)
print("wrote", out)
print(json.dumps(params, indent=2))
