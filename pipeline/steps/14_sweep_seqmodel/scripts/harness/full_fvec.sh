#!/usr/bin/env bash
# Subsample the aggregated combined ms + fvec (empirical vcfForMask) for a 2000/class
# diploSHIC model. Focal (hard_5,soft_5) fvec'd fully (2000); linked positions subsampled
# to 300 (makeTrainingSets balances to 2000/class from the 10-position pool); neutral to 3000.
# ~13000 reps. Run detached:  setsid bash full_fvec.sh > <log> 2>&1 </dev/null &
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
H=$D/scripts/harness; C=$H/combined
FR=$D/results/fullrun; FV=$FR/fvec; SUB=$FR/subms; LOGD=$FR/fveclogs
mkdir -p "$FV" "$SUB" "$LOGD"
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
PY=/home/ssmall/miniforge3/envs/slim_sim/bin/python
VCF=$D/results/task17_maskval/cleaned.1.vcf.gz
S2P=$D/results/task17_maskval/s2p.tsv
MASK=$H/talapas/mask/mask.fa   # chr1-only (106MB) — avoids 40x3.3G-mask RAM blowup
NHAP=700
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export SUB FV LOGD DE VCF S2P MASK

names="neutral hard_5 soft_5"
for i in 0 1 2 3 4 6 7 8 9 10; do names="$names hard_$i soft_$i"; done

echo "FULLFVEC_START $(date -u +%FT%TZ)"
# Right-sized $SUB/<name>.sub.msOut.gz are produced upstream by run_all.sh's aggregate
# (combine only the needed tasks). fvec each (parallel, atomic mkdir lock; fveclogs
# OUT of the fvec dir so makeTrainingSets' <prefix>_<i>.* glob can't eat them).
do_one() {
  local name=$1 ms=$SUB/$1.sub.msOut.gz out=$FV/$1.fvec lock=$FV/$1.lock
  [ -s "$ms" ]  || { echo "SKIP $name (no ms)"; return; }
  [ -s "$out" ] && { echo "SKIP $name (done)"; return; }
  mkdir "$lock" 2>/dev/null || { echo "SKIP $name (locked)"; return; }
  if "$DE/diploSHIC" fvecSim diploid "$ms" "$out.tmp" --totalPhysLen 1100000 --numSubWins 11 \
        --maskFileName "$MASK" --chrArmsForMasking 1 --vcfForMaskFileName "$VCF" \
        --popForMask illex --sampleToPopFileName "$S2P" \
        --unmaskedGenoFracCutoff 0.5 --unmaskedFracCutoff 0.25 > "$LOGD/$name.fveclog" 2>&1; then
    mv "$out.tmp" "$out"; echo "DONE $name ($(($(wc -l <"$out")-1)) rows)"
  else rm -f "$out.tmp"; echo "FAIL $name (see $LOGD/$name.fveclog)"; fi
  rmdir "$lock" 2>/dev/null
}
export -f do_one
rm -rf "$FV"/*.lock 2>/dev/null
printf '%s\n' $names | xargs -P 40 -I{} bash -c 'do_one "$@"' _ {}
rm -rf "$FV"/*.lock 2>/dev/null
echo "FULLFVEC_DONE $(date -u +%FT%TZ)"
