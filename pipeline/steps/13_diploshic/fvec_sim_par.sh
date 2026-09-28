#!/usr/bin/env bash
# Scaled rebuild of diploSHIC training feature vectors from the 44k balanced sims.
#   [1] split_sweep_hardsoft.py  -> fvec/msc/{hard,soft}_$s.msc + neutral.msc  (TRUE hard/soft labels from raw -f)
#   [2] fvecSim on all 23 msc files IN PARALLEL (P12, BLAS-thread-capped)      -> fvec/simstats/*.fvec
#   [3] makeTrainingSets (module called DIRECTLY -- CLI wrapper rejects comma lists) -> fvec/train
#   [4] split_traintest.py (80/20)
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh
DENV=/home/ssmall/miniforge3/envs/diploshic_env
PYBIN=$DENV/bin/python
AUT_CSV="${AUT// /,}"
PAR=${1:-12}
mkdir -p $WD/fvec/msc $WD/fvec/simstats $WD/fvec/train $WD/fvec/test $WD/logs

echo "[$(date)] [1] split hard/soft by TRUE label (raw discoal -f)"
$PYBIN $WD/split_sweep_hardsoft.py 2>&1 | tee $WD/logs/split_hardsoft.log
[ ! -s $WD/fvec/msc/neutral.msc ] && { echo "SPLIT FAIL: no neutral.msc"; exit 1; }

echo "[$(date)] [2] parallel fvecSim (P$PAR)"
run_fvec(){
  local name=$1
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  $DENV/bin/diploSHIC fvecSim diploid $WD/fvec/msc/${name}.msc $WD/fvec/simstats/${name}.fvec \
    --totalPhysLen 44000 --numSubWins 11 --maskFileName $WD/mask.fa \
    --chrArmsForMasking $AUT_CSV --unmaskedFracCutoff 0.25 \
    > $WD/logs/fvecsim_${name}.log 2>&1 && echo "  ok $name ($(($(wc -l < $WD/fvec/simstats/${name}.fvec)-1)) reps)" || echo "  FAIL $name"
}
export -f run_fvec; export WD DENV AUT_CSV
( echo neutral; for s in $(seq 0 10); do echo hard_$s; echo soft_$s; done ) \
  | xargs -P$PAR -I{} bash -c 'run_fvec "$@"' _ {}

# guard: every expected fvec present and non-empty
for name in neutral $(for s in $(seq 0 10); do echo hard_$s soft_$s; done); do
  [ -s $WD/fvec/simstats/${name}.fvec ] || { echo "FVEC MISSING: $name"; exit 1; }
done

echo "[$(date)] [3] makeTrainingSets (direct module call)"
MTS=$($PYBIN -c "import diploshic, os; print(os.path.join(os.path.dirname(diploshic.__file__), 'makeTrainingSets.py'))")
$PYBIN "$MTS" \
  $WD/fvec/simstats/neutral.fvec \
  $WD/fvec/simstats/soft \
  $WD/fvec/simstats/hard \
  5 \
  0,1,2,3,4,6,7,8,9,10 \
  $WD/fvec/train 2>&1 | tee $WD/logs/makeTrainingSets.log

echo "[$(date)] [4] 80/20 train/test split"
$PYBIN $WD/split_traintest.py 2>&1 | tee $WD/logs/split_traintest.log

echo "[$(date)] FVEC_PIPE_DONE"
echo "=== train ==="; wc -l $WD/fvec/train/*.fvec 2>/dev/null
echo "=== test ===";  wc -l $WD/fvec/test/*.fvec 2>/dev/null
