#!/usr/bin/env bash
# 3-GPU PARALLEL ReLERNN with THREAD CAPS.
# Fix: TF/BLAS default to 256-wide thread pools on this 256-core box -> futex thrash,
# GPU idles at 0%. Cap all thread pools to 8 so the input pipeline actually feeds the GPU.
# Fan the 4 tags across the 3 idle A100s: gpu0=male_auto, gpu1=female_auto, gpu2=Z_male->Z_female.
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
B=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin
ENV=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env
RD=$A/steps/11_relernn; DEM=$RD/demHist.stairway.txt; L=$RD/logs; mkdir -p $L

# thread caps (the whole point): stop TF/numpy/BLAS from spawning nproc-wide pools
export PYTHONNOUSERSITE=1 TF_USE_LEGACY_KERAS=1 TF_FORCE_GPU_ALLOW_GROWTH=true \
       XLA_FLAGS=--xla_gpu_cuda_data_dir=$ENV \
       OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
       NUMEXPR_NUM_THREADS=8 TF_NUM_INTRAOP_THREADS=8 TF_NUM_INTEROP_THREADS=2

pipe(){ local tag=$1 vcf=$2 gbed=$3 mask=$4 gpu=$5; local d=$RD/run_$tag
  export CUDA_VISIBLE_DEVICES=$gpu
  find $d/proj -name '*BSCORRECTED.txt' 2>/dev/null | grep -q . && { echo "[$tag] already done"; return 0; }
  # SIMULATE once (only if sims absent) -- current run has these cached, so skipped
  if [ ! -f $d/proj/train/info.p ] || [ ! -f $d/proj/vali/info.p ] || [ ! -f $d/proj/test/info.p ]; then
    rm -rf $d/proj; echo "[$(date)] [$tag] SIMULATE (gpu$gpu)"
    $B/ReLERNN_SIMULATE -v $vcf -g $gbed -m $mask -d $d/proj -n $DEM -u 3e-9 -l 1 -t 12 \
      --unphased --maskThresh 0.9 --nTrain 12000 --nVali 1000 --nTest 1000 > $L/$tag.sim.log 2>&1
    { [ -f $d/proj/train/info.p ] && [ -f $d/proj/vali/info.p ] && [ -f $d/proj/test/info.p ]; } || { echo "[$tag] SIM FAIL"; grep -iE 'error|zerodiv|both haploid' $L/$tag.sim.log|tail -2; return 1; }
  else echo "[$tag] reuse existing sims (gpu$gpu)"; fi
  # TRAIN with retry (do NOT rm proj -> keep sims)
  local ok=0
  for a in 1 2 3 4; do
    echo "[$(date)] [$tag] TRAIN attempt $a (gpu$gpu)"
    $B/ReLERNN_TRAIN -d $d/proj -t 8 --nEpochs 50 --nValSteps 20 > $L/$tag.train.log 2>&1 && { ok=1; break; }
    echo "[$tag] TRAIN died (attempt $a), retrying"; sleep 20
  done
  [ $ok -eq 0 ] && { echo "[$tag] TRAIN FAILED x4"; return 1; }
  echo "[$(date)] [$tag] PREDICT"; for a in 1 2 3; do $B/ReLERNN_PREDICT -v $vcf -d $d/proj --unphased > $L/$tag.pred.log 2>&1 && break; sleep 15; done
  echo "[$(date)] [$tag] BSCORRECT"; for a in 1 2 3; do $B/ReLERNN_BSCORRECT -d $d/proj -t 8 --nSlice 20 --nReps 20 > $L/$tag.bs.log 2>&1 && break; sleep 15; done
  local out=$(find $d/proj -name '*BSCORRECTED.txt'|head -1)
  [ -n "$out" ] && awk 'NR==1{for(i=1;i<=NF;i++)if($i=="recombRate")r=i} NR>1{s+=$r;n++} END{if(n)printf "[%s] mean r=%.3e/bp = %.4f cM/Mb (n=%d)\n","'$tag'",s/n,s/n*1e8,n}' "$out"
}

# per-GPU workers (each its own subshell -> isolated CUDA_VISIBLE_DEVICES). Z is 1 chrom (fast) so both Z tags share gpu2.
gpu0(){ pipe male_auto   $RD/male/male.kept.vcf     $RD/genome.bed        $RD/inaccessible.bed 0; }
gpu1(){ pipe female_auto $RD/female/female.kept.vcf $RD/genome.bed        $RD/inaccessible.bed 1; }
gpu2(){ pipe Z_male      $RD/Z_male/Z_male.kept.vcf $RD/Z_male/genome.bed $RD/Z_male/inacc.bed 2
        pipe Z_female    $RD/Z_female/Z_female.kept.vcf $RD/Z_female/genome.bed $RD/Z_female/inacc.bed 2; }

echo "[$(date)] PARALLEL 3-GPU ReLERNN start"
gpu0 > $L/gpu0.worker.log 2>&1 &
gpu1 > $L/gpu1.worker.log 2>&1 &
gpu2 > $L/gpu2.worker.log 2>&1 &
wait
echo "[$(date)] ALL RUNS DONE"
