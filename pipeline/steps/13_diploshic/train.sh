#!/usr/bin/env bash
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh
export CUDA_DEVICE_ORDER=PCI_BUS_ID
# TF thread caps (shared 256-core box — avoid futex thrash, see relernn-thread-cap memory)
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export TF_NUM_INTRAOP_THREADS=8 TF_NUM_INTEROP_THREADS=2
# This net is tiny (132-feature dense CNN, ~2k examples) and trains in seconds.
# The diploshic_env TF cannot init the CUDA driver in this invocation (falls back to
# CPU regardless), so pin to CPU explicitly to avoid a partial GPU grab on the shared box.
export CUDA_VISIBLE_DEVICES=""
$DSHIC train $WD/fvec/train/ $WD/fvec/test/ $WD/model \
  --epochs 100 --confusionFile $WD/confusion_train.png \
  2>&1 | tee $WD/train.log
# diploSHIC 1.x writes model.weights.h5 (not .hdf5)
ls $WD/model.json $WD/model.weights.h5
