#!/usr/bin/env bash
set -uo pipefail
export PYTHONNOUSERSITE=1 CUDA_VISIBLE_DEVICES=0 TF_USE_LEGACY_KERAS=1
export XLA_FLAGS=--xla_gpu_cuda_data_dir=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env
cd /sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn/pilot_chr45_male
echo "[$(date)] TRAIN"; /home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin/ReLERNN_TRAIN -d proj -t 12 --nEpochs 50 --nValSteps 10 > train4.log 2>&1 && echo TRAIN_OK || echo TRAIN_FAIL
echo "[$(date)] PREDICT"; /home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin/ReLERNN_PREDICT -v chr45_male.kept.vcf -d proj -u 3e-9 -l 1 > predict4.log 2>&1 && echo PRED_OK || echo PRED_FAIL
echo "[$(date)] BSCORRECT"; /home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin/ReLERNN_BSCORRECT -d proj -t 12 --nSlice 20 --nReps 20 > bscorrect4.log 2>&1 && echo BS_OK || echo BS_FAIL
echo "[$(date)] DONE"
