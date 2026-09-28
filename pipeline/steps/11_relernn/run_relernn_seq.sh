#!/usr/bin/env bash
# SEQUENTIAL ReLERNN on GPU0 ONLY, with THREAD CAPS.
# Rationale: shared GPU box -- other users grab GPU1/GPU2, and each ReLERNN model needs ~74GB
# (nearly a full A100), so parallel-across-cards OOMs ("Dst tensor not initialized") when a
# neighbor card isn't actually free. We reliably HOLD GPU0 (74GB occupancy keeps others out).
# One training at a time also gives that job the full CPU batch-prep pipeline -> faster steps.
# Thread caps stop TF/BLAS 256-wide pool futex-thrash. Reuses all cached sims.
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
B=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin
ENV=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env
RD=$A/steps/11_relernn; DEM=$RD/demHist.stairway.txt; L=$RD/logs; mkdir -p $L
export PYTHONNOUSERSITE=1 TF_USE_LEGACY_KERAS=1 TF_FORCE_GPU_ALLOW_GROWTH=true \
       XLA_FLAGS=--xla_gpu_cuda_data_dir=$ENV CUDA_VISIBLE_DEVICES=0 \
       OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
       NUMEXPR_NUM_THREADS=8 TF_NUM_INTRAOP_THREADS=8 TF_NUM_INTEROP_THREADS=2

pipe(){ local tag=$1 vcf=$2 gbed=$3 mask=$4; local d=$RD/run_$tag
  find $d/proj -name '*BSCORRECTED.txt' 2>/dev/null | grep -q . && { echo "[$(date)] [$tag] already done"; return 0; }
  if [ ! -f $d/proj/train/info.p ] || [ ! -f $d/proj/vali/info.p ] || [ ! -f $d/proj/test/info.p ]; then
    rm -rf $d/proj; echo "[$(date)] [$tag] SIMULATE"
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    $B/ReLERNN_SIMULATE -v $vcf -g $gbed -m $mask -d $d/proj -n $DEM -u 3e-9 -l 1 -t 12 \
      --unphased --maskThresh 0.9 --nTrain 12000 --nVali 1000 --nTest 1000 > $L/$tag.sim.log 2>&1
    { [ -f $d/proj/train/info.p ] && [ -f $d/proj/vali/info.p ] && [ -f $d/proj/test/info.p ]; } || { echo "[$tag] SIM FAIL"; grep -iE 'error|zerodiv|both haploid' $L/$tag.sim.log|tail -2; return 1; }
  else echo "[$(date)] [$tag] reuse existing sims"; fi
  local ok=0
  for a in 1 2 3 4; do
    echo "[$(date)] [$tag] TRAIN attempt $a (gpu0)"
    $B/ReLERNN_TRAIN -d $d/proj -t 8 --nEpochs 50 --nValSteps 20 > $L/$tag.train.log 2>&1 && { ok=1; break; }
    echo "[$(date)] [$tag] TRAIN died (attempt $a), retrying"; sleep 30
  done
  [ $ok -eq 0 ] && { echo "[$tag] TRAIN FAILED x4"; return 1; }
  echo "[$(date)] [$tag] PREDICT"; for a in 1 2 3; do $B/ReLERNN_PREDICT -v $vcf -d $d/proj --unphased --batchSizeOverride 50 > $L/$tag.pred.log 2>&1 && break; sleep 15; done
  echo "[$(date)] [$tag] BSCORRECT"; for a in 1 2 3; do $B/ReLERNN_BSCORRECT -d $d/proj -t 8 --nSlice 20 --nReps 20 > $L/$tag.bs.log 2>&1 && break; sleep 15; done
  local out=$(find $d/proj -name '*BSCORRECTED.txt'|head -1)
  [ -n "$out" ] && awk 'NR==1{for(i=1;i<=NF;i++)if($i=="recombRate")r=i} NR>1{s+=$r;n++} END{if(n)printf "[%s] mean r=%.3e/bp = %.4f cM/Mb (n=%d)\n","'$tag'",s/n,s/n*1e8,n}' "$out"
}
echo "[$(date)] SEQUENTIAL(GPU0) ReLERNN start"
pipe male_auto   $RD/male/male.kept.vcf     $RD/genome.bed         $RD/inaccessible.bed
pipe female_auto $RD/female/female.kept.vcf $RD/genome.bed         $RD/inaccessible.bed
pipe Z_male      $RD/Z_male/Z_male.kept.vcf $RD/Z_male/genome.bed  $RD/Z_male/inacc.bed
pipe Z_female    $RD/Z_female/Z_female.kept.vcf $RD/Z_female/genome.bed $RD/Z_female/inacc.bed
echo "[$(date)] ALL SEQUENTIAL RUNS DONE"
