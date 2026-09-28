#!/usr/bin/env bash
# Generate additional discoal sweep sims, UNIFORM across all 11 subwindows, to
# reach ~4000 reps/window (target ~2000 balanced examples/class). makeTrainingSets
# draws the hard/soft classes from the CENTER window only and balances all 5
# classes to the min, so uniform allocation needs ~2x the per-class target per
# window. Reuses run_one.py (S.sweep: ~50/50 hard/soft per rep). The TRUE
# hard/soft label is recovered later from the raw discoal command line ('-f' =>
# soft) by split_sweep_hardsoft.py -- NOT from any RNG replay.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic
PYBIN=/home/ssmall/miniforge3/envs/diploshic_env/bin/python
export DISCOAL=/home/ssmall/programs/discoal_oobfix/discoal
mkdir -p "$D/sims/raw/sweep"

REP_START=${1:-201}
REP_END=${2:-4000}
PAR=${3:-48}

gen(){
  local sub=$1 rep=$2
  local out="$D/sims/raw/sweep/sweep_s${sub}_r${rep}.ms"
  [ -s "$out" ] && return 0                 # resumable: skip existing non-empty
  "$PYBIN" "$D/run_one.py" sweep "$sub" "$rep" "$out"
}
export -f gen; export D PYBIN DISCOAL

echo "simulate_more: sweep reps r${REP_START}..r${REP_END} across subwin 0..10, P${PAR}"
( for s in $(seq 0 10); do for r in $(seq "$REP_START" "$REP_END"); do echo "$s $r"; done; done ) \
  | xargs -P"$PAR" -n2 bash -c 'gen "$@"' _

echo "DONE simulate_more: sweep files now $(ls "$D"/sims/raw/sweep/*.ms | wc -l)"
