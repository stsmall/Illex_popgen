#!/usr/bin/env bash
# Parallel diploSHIC fvecSim for the pilot, with the empirical vcfForMask workflow.
# Concurrency-safe skip via an atomic mkdir lock (diploSHIC creates its output only
# late, so a .tmp-existence guard does NOT prevent a second launch from re-running an
# in-progress file). Writes to .tmp then renames. Run as ONE invocation.
#   bash fvec_pilot.sh [CONCURRENCY]
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
OUT=$D/results/pilot/train_ms
FV=$D/results/pilot/fvec; mkdir -p "$FV"
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
VCF=$D/results/task17_maskval/cleaned.1.vcf.gz
S2P=$D/results/task17_maskval/s2p.tsv
MASK=$D/../13_diploshic/mask.fa
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export OUT FV DE VCF S2P MASK
CONC=${1:-12}

do_one() {
  local name=$1 ms=$OUT/$1.msOut.gz out=$FV/$1.fvec lock=$FV/$1.lock
  [ -f "$ms" ]  || { echo "SKIP $name (no ms yet)"; return; }
  [ -s "$out" ] && { echo "SKIP $name (fvec done)"; return; }
  mkdir "$lock" 2>/dev/null || { echo "SKIP $name (locked/in progress)"; return; }
  if "$DE/diploSHIC" fvecSim diploid "$ms" "$out.tmp" \
        --totalPhysLen 1100000 --numSubWins 11 \
        --maskFileName "$MASK" --chrArmsForMasking 1 --vcfForMaskFileName "$VCF" \
        --popForMask illex --sampleToPopFileName "$S2P" \
        --unmaskedGenoFracCutoff 0.5 --unmaskedFracCutoff 0.25 \
        > "$FV/$name.fveclog" 2>&1; then
    mv "$out.tmp" "$out"; echo "DONE $name ($(($(wc -l <"$out")-1)) rows)"
  else
    rm -f "$out.tmp"; echo "FAIL $name (see $name.fveclog)"
  fi
  rmdir "$lock" 2>/dev/null
}
export -f do_one

# worklist: neutral + focal/linked positions we actually have (hard_9/10 lost)
names="neutral"
for i in 0 1 2 3 4 5 6 7 8; do names="$names hard_$i soft_$i"; done
rm -rf "$FV"/*.lock 2>/dev/null      # clear stale locks from any killed run
echo "FVEC_PILOT_START $(date -u +%H:%M:%S) conc=$CONC"
printf '%s\n' $names | xargs -P "$CONC" -I{} bash -c 'do_one "$@"' _ {}
echo "FVEC_PILOT_DONE $(date -u +%H:%M:%S)"
