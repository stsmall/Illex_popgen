#!/usr/bin/env bash
# Predict the full-SFS retrained model (illexModel_vcf) on the 43 pre-computed empirical full-SFS fvecs.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
MODEL=$D/results/vcfretrain_fullsfs/illexModel_vcf
FVD=$D/results/empirical_scan_fullsfs/fvec
OUT=$D/results/empirical_scan_fullsfs; PRED=$OUT/preds; mkdir -p "$PRED" "$OUT/predlogs"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 CUDA_VISIBLE_DEVICES="" TF_CPP_MIN_LOG_LEVEL=3
export DE MODEL FVD PRED OUT
PROG=$OUT/predict.progress; : > "$PROG"
do_pred(){
  local f=$1 c; c=$(basename "$f" .fvec); local out=$PRED/$c.preds
  [ -s "$out" ] && { echo "skip $c"; return; }
  "$DE/diploSHIC" predict "$MODEL.json" "$MODEL.weights.h5" "$f" "$out" --numSubWins 11 \
     > "$OUT/predlogs/$c.log" 2>&1 \
   && echo "$(date -u +%FT%TZ) DONE $c ($(($(wc -l <"$out")-1)) win)" >> "$PROG" \
   || echo "$(date -u +%FT%TZ) FAIL $c" >> "$PROG"
}
export -f do_pred
echo "$(date -u +%FT%TZ) PREDICT_START ($(ls $FVD/*.fvec|wc -l) fvecs)" >> "$PROG"
ls "$FVD"/*.fvec | xargs -P 6 -I{} bash -c 'do_pred "$@"' _ {}
np=$(ls "$PRED"/*.preds 2>/dev/null|wc -l)
echo "$(date -u +%FT%TZ) PREDICT_DONE ($np/43)" >> "$PROG"
# combine
first=$(ls "$PRED"/*.preds 2>/dev/null|head -1)
{ head -1 "$first"; for f in "$PRED"/*.preds; do tail -n +2 "$f"; done; } > "$OUT/genome.preds"
echo "$(date -u +%FT%TZ) COMBINED genome.preds ($(($(wc -l <"$OUT/genome.preds")-1)) windows)" >> "$PROG"
