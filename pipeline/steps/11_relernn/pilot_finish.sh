#!/usr/bin/env bash
# Wait for pilot SIMULATE to finish, then TRAIN -> PREDICT -> BSCORRECT (chr45 male, equilibrium pilot).
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
P=$A/steps/11_relernn/pilot_chr45_male
B=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin
export PYTHONNOUSERSITE=1 CUDA_VISIBLE_DEVICES=0
# wait until SIMULATE done: train/ has .npz sims and no ReLERNN_SIMULATE proc
while pgrep -f "ReLERNN_SIMULATE.*pilot_chr45" >/dev/null; do sleep 15; done
echo "[$(date)] SIMULATE done: train=$(ls $P/proj/train 2>/dev/null|wc -l) vali=$(ls $P/proj/vali 2>/dev/null|wc -l)"
echo "[$(date)] TRAIN"
$B/ReLERNN_TRAIN -d $P/proj -t 12 --nEpochs 50 --nValSteps 10 > $P/train.log 2>&1
echo "[$(date)] PREDICT"
$B/ReLERNN_PREDICT -v $P/chr45_male.vcf -d $P/proj -u 3e-9 -l 1 > $P/predict.log 2>&1
echo "[$(date)] BSCORRECT"
$B/ReLERNN_BSCORRECT -d $P/proj -t 12 --nSlice 20 --nReps 20 > $P/bscorrect.log 2>&1
echo "[$(date)] PILOT DONE"
ls -la $P/proj/*.txt $P/proj/*PREDICT* 2>/dev/null | head
