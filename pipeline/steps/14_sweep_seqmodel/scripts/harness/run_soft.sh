#!/usr/bin/env bash
# Soft-only pilot campaign (hard_0-8 + neutral salvaged from the interrupted run).
# Detached so it survives CLI idle:  setsid bash run_soft.sh > <log> 2>&1 </dev/null &
set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel
RC=$D/scripts/harness/run_campaign.py
OUT=$D/results/pilot/train_ms
PY=/home/ssmall/miniforge3/envs/slim_sim/bin/python
echo "SOFT_START $(date -u +%H:%M:%S)"
$PY $RC --klass soft --n 40 --concurrency 16 --outdir "$OUT" --seed0 3000
echo "SOFT_DONE $(date -u +%H:%M:%S)"
ls -la "$OUT"/soft_*.msOut.gz
echo "PILOT_GEN_DONE $(date -u +%H:%M:%S)"
