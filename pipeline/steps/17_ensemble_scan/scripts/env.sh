#!/usr/bin/env bash
# Shared paths/config for the illex ensemble sweep scan.
B=/sietch_colab/data_share/illex/popgen_data
E=$B/analysis/steps/17_ensemble_scan
S=$E/scripts; R=$E/results; L=$E/logs
BCF=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
export BCFTOOLS_PLUGINS=/home/ssmall/miniforge3/envs/bioinfo-buddy/libexec/bcftools
PY=/home/ssmall/miniforge3/envs/slim_sim/bin/python
REF=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa
ANC=/sietch_colab/data_share/illex/assembly/current_assembly_version/ancestral_ref.fa
SF2=/home/ssmall/programs/SweepFinder2/SweepFinder2
RAISD=/home/ssmall/programs/RAiSD/bin/release/RAiSD          # rebuilt with zlib (reads .vcf.gz)
BVAL=/home/ssmall/miniforge3/envs/bvalcalc/bin/bvalcalc
BVALDIR=$B/analysis/steps/16_bgspy/bvalcalc
# SOURCE callset = baker_2025 (FULL SFS: singletons + high-freq-derived; NOT the MAF-filtered
# 13_diploshic/combined, which is gutted for SFS-based tests).
VCFDIR=$B/seq_data/baker_2025/vcfs
N0=200                             # SweepFinder2 constant projected sample size (2*100 diploids)
GRID_BP=100000                     # SF2 grid spacing (bp); RAiSD uses a finer 10kb grid
CHROMS="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"
