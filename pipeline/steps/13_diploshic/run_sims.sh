#!/usr/bin/env bash
# Race-safe: parallel phase writes only .ms files; manifest is built afterward
# by listing produced files and parsing kind/subwin from filenames.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh
NREP=${1:-2000}; mkdir -p $WD/sims/raw/neutral $WD/sims/raw/sweep

gen(){
  local kind=$1 sub=$2 rep=$3
  local out=$WD/sims/raw/$kind/${kind}_s${sub}_r${rep}.ms
  $PY $WD/run_one.py "$kind" "$sub" "$rep" "$out"
}
export -f gen; export WD PY DISCOAL

# Run parallel phase — each job writes only its own .ms file (no shared-manifest writes)
( for r in $(seq 1 $NREP); do echo "neutral -1 $r"; done
  for s in $(seq 0 10); do for r in $(seq 1 $((NREP/11))); do echo "sweep $s $r"; done; done ) \
  | xargs -P16 -I{} bash -c 'gen $0' "{}"

# Build manifest after all parallel jobs complete (race-free)
echo -e "path\tkind\tsweep_subwin" > $WD/sims/manifest.tsv
# Neutral: subwin encoded as -1 in filename as "s-1"
for f in $WD/sims/raw/neutral/neutral_s-1_r*.ms; do
    [ -f "$f" ] && echo -e "$f\tneutral\t-1"
done >> $WD/sims/manifest.tsv
# Sweeps: subwin is 0..10
for s in $(seq 0 10); do
    for f in $WD/sims/raw/sweep/sweep_s${s}_r*.ms; do
        [ -f "$f" ] && echo -e "$f\tsweep\t$s"
    done
done >> $WD/sims/manifest.tsv

echo "sims done: $(( $(wc -l < $WD/sims/manifest.tsv) - 1 )) rows (excluding header)"
