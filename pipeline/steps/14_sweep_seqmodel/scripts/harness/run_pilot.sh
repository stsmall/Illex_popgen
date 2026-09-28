#!/usr/bin/env bash
# Pilot training-sim generation, run DETACHED (setsid) so it survives CLI session
# idle. Sequential campaigns; progress + sentinels to the log so it can be polled.
#   setsid bash run_pilot.sh > <log> 2>&1 < /dev/null &
set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
RC=$D/scripts/harness/run_campaign.py
OUT=$D/results/pilot/train_ms
PY=/home/ssmall/miniforge3/envs/slim_sim/bin/python

echo "PILOT_START $(date -u +%H:%M:%S)"
$PY $RC --klass neutral --n 6  --concurrency 8  --outdir "$OUT" --seed0 1000
echo "NEUTRAL_DONE $(date -u +%H:%M:%S)"
$PY $RC --klass hard    --n 60 --concurrency 16 --outdir "$OUT" --seed0 2000
echo "HARD_DONE $(date -u +%H:%M:%S)"
$PY $RC --klass soft    --n 60 --concurrency 16 --outdir "$OUT" --seed0 3000
echo "SOFT_DONE $(date -u +%H:%M:%S)"
echo "=== combined training ms ==="; ls -la "$OUT"/*.msOut.gz
echo "PILOT_GEN_DONE $(date -u +%H:%M:%S)"
