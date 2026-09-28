#!/usr/bin/env bash
# Genome-wide male PREDICT using the backported chunked ReLERNN (v2.0 --batchSizeOverride logic).
# Memory-safe (per-batch, GT cached per-chrom). Uses salvaged trained weights. Then BSCORRECT.
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
B=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin
ENV=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env
RD=$A/steps/11_relernn; d=$RD/run_male_auto; proj=$d/proj; L=$RD/logs
export PYTHONNOUSERSITE=1 TF_USE_LEGACY_KERAS=1 TF_FORCE_GPU_ALLOW_GROWTH=true \
       XLA_FLAGS=--xla_gpu_cuda_data_dir=$ENV CUDA_VISIBLE_DEVICES=0 \
       OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
       NUMEXPR_NUM_THREADS=8 TF_NUM_INTRAOP_THREADS=8 TF_NUM_INTEROP_THREADS=2
# restore salvaged weights + full 43-chrom windowSizes
/bin/cp -pf $proj/networks/weights.attempt1_ep50.h5 $proj/networks/weights.h5
/bin/cp -pf $proj/networks/model.attempt1.json $proj/networks/model.json
/bin/cp -pf $proj/networks/windowSizes.txt.full $proj/networks/windowSizes.txt
echo "[$(date)] male genome-wide PREDICT start ($(wc -l <$proj/networks/windowSizes.txt) chroms, bso=50)"
for a in 1 2 3; do
  $B/ReLERNN_PREDICT -v $RD/male/male.kept.vcf -d $proj --unphased --batchSizeOverride 50 > $L/male_auto.pred.log 2>&1 && break
  echo "[$(date)] PREDICT attempt $a exited $? ; retry"; sleep 20
done
if [ -f $proj/male.kept.PREDICT.txt ]; then
  echo "[$(date)] male PREDICT done: $(wc -l <$proj/male.kept.PREDICT.txt) windows"
  awk 'NR>1{s+=$5;n++}END{if(n)printf "[male PREDICT] mean=%.4f cM/Mb n=%d\n",s/n*1e8,n}' $proj/male.kept.PREDICT.txt
  echo "[$(date)] male BSCORRECT start"
  for a in 1 2 3; do $B/ReLERNN_BSCORRECT -d $proj -t 8 --nSlice 20 --nReps 20 > $L/male_auto.bs.log 2>&1 && break; sleep 20; done
  out=$(find $proj -name '*BSCORRECTED.txt'|head -1)
  echo "[$(date)] male BSCORRECTED: ${out:-MISSING}"
  [ -n "$out" ] && awk 'NR==1{for(i=1;i<=NF;i++)if($i=="recombRate")r=i}NR>1{s+=$r;n++}END{if(n)printf "[male BSCORRECTED] mean=%.4f cM/Mb n=%d\n",s/n*1e8,n}' "$out"
else
  echo "[$(date)] male PREDICT FAILED (no PREDICT.txt) - check $L/male_auto.pred.log"
fi
echo "[$(date)] male_predict_gw DONE"
