#!/usr/bin/env bash
# SALVAGE Z_male: attempt-1 TRAIN completed all 50 epochs and saved the model
# (networks/weights.h5, model.json @16:49), then OOM'd in the OPTIONAL post-training
# test-set eval (model.predict on full test tensor; Z_male's 658 haps > Z_female's 560).
# The driver blindly started a full retrain (attempt 2) which would hit the SAME
# deterministic OOM after ~28h. Killed it; model is intact. Run PREDICT (chunked,
# won't OOM) + BSCORRECT directly on the saved model. Mirrors run_one_tag.sh 6-20,37-42.
set -uo pipefail
tag=Z_male; gpu=0
A=/sietch_colab/data_share/illex/popgen_data/analysis
B=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin
ENV=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env
RD=$A/steps/11_relernn; L=$RD/logs
d=$RD/run_$tag; vcf=$RD/$tag/$tag.kept.vcf
gpuuuid=$(nvidia-smi --query-gpu=uuid --format=csv,noheader -i $gpu 2>/dev/null)
echo "[$(date)] [$tag] SALVAGE predict on GPU$gpu ($gpuuuid)"
export PYTHONNOUSERSITE=1 TF_USE_LEGACY_KERAS=1 TF_FORCE_GPU_ALLOW_GROWTH=true \
       XLA_FLAGS=--xla_gpu_cuda_data_dir=$ENV CUDA_DEVICE_ORDER=PCI_BUS_ID \
       CUDA_VISIBLE_DEVICES=${gpuuuid:-$gpu} \
       OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
       NUMEXPR_NUM_THREADS=8 TF_NUM_INTRAOP_THREADS=8 TF_NUM_INTEROP_THREADS=2
echo "[$(date)] [$tag] PREDICT"; for a in 1 2 3; do $B/ReLERNN_PREDICT -v $vcf -d $d/proj --unphased --batchSizeOverride 50 --gpuID $gpu > $L/$tag.pred.log 2>&1 && break; echo "[$tag] predict retry $a"; sleep 15; done
echo "[$(date)] [$tag] BSCORRECT"; for a in 1 2 3; do $B/ReLERNN_BSCORRECT -d $d/proj -t 8 --nSlice 20 --nReps 20 --gpuID $gpu > $L/$tag.bs.log 2>&1 && break; echo "[$tag] bscorrect retry $a"; sleep 15; done
out=$(find $d/proj -name '*BSCORRECTED.txt'|head -1)
echo "[$(date)] [$tag] BSCORRECTED: ${out:-MISSING}"
[ -n "$out" ] && awk -v tag=$tag 'NR==1{for(i=1;i<=NF;i++)if($i=="recombRate")r=i}NR>1{s+=$r;n++}END{if(n)printf "[%s] mean=%.4f cM/Mb (n=%d)\n",tag,s/n*1e8,n}' "$out"
echo "[$(date)] [$tag] SALVAGE DONE"
