#!/usr/bin/env bash
export ANA=/sietch_colab/data_share/illex/popgen_data/analysis
export WD=$ANA/steps/13_diploshic
export DENV=/home/ssmall/miniforge3/envs/diploshic_env
export DISCOAL=/home/ssmall/programs/discoal_oobfix/discoal
export DSHIC="mamba run -p $DENV diploSHIC"
export PY="mamba run -p $DENV python"
export BB=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
export REF=/sietch_colab/data_share/illex/andy_links/Illex_F24.primary.clean.fa
export AUT="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"
