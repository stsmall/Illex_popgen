#!/usr/bin/env bash
# Run ONE ReLERNN tag end-to-end on a SPECIFIED GPU (parallel Z runs on GPU1/GPU2).
# Usage: run_one_tag.sh <tag> <gpuid>   (tag e.g. Z_male, Z_female). Uses anti-leak chunked predict.
set -uo pipefail
tag=$1; gpu=$2
A=/sietch_colab/data_share/illex/popgen_data/analysis
B=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin
ENV=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env
RD=$A/steps/11_relernn; L=$RD/logs; DEM=$RD/demHist.stairway.txt
d=$RD/run_$tag; vcf=$RD/$tag/$tag.kept.vcf; gbed=$RD/$tag/genome.bed; mask=$RD/$tag/inacc.bed
# Pin to the target card by UUID (ordering-immune). Default CUDA "fastest-first" order does NOT
# match nvidia-smi's PCI index, so a bare CVD=<index> can land on the wrong physical GPU (e.g. the
# busy GPU0). PCI_BUS_ID makes CUDA indices == nvidia-smi indices; UUID removes all ambiguity.
gpuuuid=$(nvidia-smi --query-gpu=uuid --format=csv,noheader -i $gpu 2>/dev/null)
echo "[$(date)] [$tag] gpu index $gpu -> $gpuuuid"
export PYTHONNOUSERSITE=1 TF_USE_LEGACY_KERAS=1 TF_FORCE_GPU_ALLOW_GROWTH=true \
       XLA_FLAGS=--xla_gpu_cuda_data_dir=$ENV CUDA_DEVICE_ORDER=PCI_BUS_ID \
       CUDA_VISIBLE_DEVICES=${gpuuuid:-$gpu} \
       OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
       NUMEXPR_NUM_THREADS=8 TF_NUM_INTRAOP_THREADS=8 TF_NUM_INTEROP_THREADS=2
echo "[$(date)] [$tag] START on GPU$gpu"
find $d/proj -name '*BSCORRECTED.txt' 2>/dev/null | grep -q . && { echo "[$tag] already done"; exit 0; }
# SIMULATE if sims incomplete (Z_female); OMP=1 to avoid fork-after-threads deadlock
if [ ! -f $d/proj/train/info.p ] || [ ! -f $d/proj/vali/info.p ] || [ ! -f $d/proj/test/info.p ]; then
  rm -rf $d/proj; echo "[$(date)] [$tag] SIMULATE"
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  $B/ReLERNN_SIMULATE -v $vcf -g $gbed -m $mask -d $d/proj -n $DEM -u 3e-9 -l 1 -t 12 \
    --unphased ${SIM_EXTRA:-} --maskThresh 0.9 --nTrain 12000 --nVali 1000 --nTest 1000 > $L/$tag.sim.log 2>&1
  { [ -f $d/proj/train/info.p ] && [ -f $d/proj/vali/info.p ] && [ -f $d/proj/test/info.p ]; } || { echo "[$tag] SIM FAIL"; grep -iE 'error|zerodiv|haploid' $L/$tag.sim.log|tail -2; exit 1; }
else echo "[$(date)] [$tag] reuse existing sims"; fi
# wait until target GPU has >=70GB free (shared box: GPU1/2 intermittently occupied by others)
wait_gpu_free(){ for w in $(seq 1 120); do free=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits -i $gpu 2>/dev/null); [ "${free:-0}" -ge 70000 ] && return 0; echo "[$(date)] [$tag] GPU$gpu only ${free}MiB free, waiting..."; sleep 60; done; return 0; }
# TRAIN (retry x8, wait for free GPU before each attempt)
ok=0; for a in $(seq 1 8); do wait_gpu_free; echo "[$(date)] [$tag] TRAIN attempt $a (gpu$gpu)"; $B/ReLERNN_TRAIN -d $d/proj -t 8 --nEpochs 50 --nValSteps 20 --gpuID $gpu > $L/$tag.train.a$a.log 2>&1 && { ok=1; /bin/cp -f $L/$tag.train.a$a.log $L/$tag.train.log; break; }; echo "[$tag] TRAIN died $a (see $L/$tag.train.a$a.log)"; /bin/cp -f $L/$tag.train.a$a.log $L/$tag.train.log; sleep 30; done
[ $ok -eq 0 ] && { echo "[$tag] TRAIN FAILED x8"; exit 1; }
# PREDICT (chunked anti-leak) then BSCORRECT
echo "[$(date)] [$tag] PREDICT"; for a in 1 2 3; do $B/ReLERNN_PREDICT -v $vcf -d $d/proj --unphased --batchSizeOverride 50 --gpuID $gpu > $L/$tag.pred.log 2>&1 && break; sleep 15; done
echo "[$(date)] [$tag] BSCORRECT"; for a in 1 2 3; do $B/ReLERNN_BSCORRECT -d $d/proj -t 8 --nSlice 20 --nReps 20 --gpuID $gpu > $L/$tag.bs.log 2>&1 && break; sleep 15; done
out=$(find $d/proj -name '*BSCORRECTED.txt'|head -1)
echo "[$(date)] [$tag] BSCORRECTED: ${out:-MISSING}"
[ -n "$out" ] && awk -v tag=$tag 'NR==1{for(i=1;i<=NF;i++)if($i=="recombRate")r=i}NR>1{s+=$r;n++}END{if(n)printf "[%s] mean=%.4f cM/Mb (n=%d)\n",tag,s/n*1e8,n}' "$out"
echo "[$(date)] [$tag] DONE"
