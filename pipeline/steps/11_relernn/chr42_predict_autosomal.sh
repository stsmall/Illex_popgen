#!/usr/bin/env bash
# Apply the CLEAN autosomal-trained ReLERNN nets (run_{male,female}_auto/proj) to chr42, BY SEX,
# via ReLERNN_PREDICT. Identical recipe to chr2_predict_autosomal.sh (chr42 is a normal autosome that
# was excluded from the ReLERNN autosomal build, like chr2). NO --forceDiploid. GPU 0 (2 is busy).
set -uo pipefail
R=/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn
A=/sietch_colab/data_share/illex/popgen_data/analysis
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
BIN=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env/bin
ENV=/home/ssmall/miniforge3/envs/grenepipe/envs/relernn_env
GPU=0
OUT=$R/chr42_autosomal_predict; mkdir -p "$OUT"; L=$R/logs; mkdir -p "$L"
P=$OUT/progress; echo "START $(date -u +%FT%TZ)" > "$P"
gpuuuid=$(nvidia-smi --query-gpu=uuid --format=csv,noheader -i $GPU 2>/dev/null)
export PYTHONNOUSERSITE=1 TF_USE_LEGACY_KERAS=1 TF_FORCE_GPU_ALLOW_GROWTH=true \
       XLA_FLAGS=--xla_gpu_cuda_data_dir=$ENV CUDA_DEVICE_ORDER=PCI_BUS_ID \
       CUDA_VISIBLE_DEVICES=${gpuuuid:-$GPU} \
       OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8 \
       TF_NUM_INTRAOP_THREADS=8 TF_NUM_INTEROP_THREADS=2

for sex in male female; do
  keep=$R/$sex/keep.txt
  proj=$R/run_${sex}_auto/proj
  vcf=$OUT/chr42_${sex}.vcf
  echo "[$sex] building chr42 VCF from the trained keep-set ($(wc -l < "$keep") indiv)" >> "$P"
  $BB/bcftools view -S "$keep" --force-samples -m2 -M2 -v snps \
      "$A/steps/00_callset/filtered/42/variants_filt.vcf.gz" -Ou 2>/dev/null \
    | $BB/bcftools view -e 'F_MISSING>0.5' -Ou 2>/dev/null \
    | $BB/bcftools +fill-tags -Ou -- -t MAF 2>/dev/null \
    | $BB/bcftools view -e 'INFO/MAF<0.01' -Oz -o "$vcf.gz" 2>/dev/null
  $BB/bgzip -dk -f "$vcf.gz"
  echo "[$sex] chr42 VCF: $($BB/bcftools index -n "$vcf.gz" 2>/dev/null) SNPs, $($BB/bcftools query -l "$vcf.gz" | wc -l) samples" >> "$P"
  echo "[$sex] PREDICT with autosomal net ($proj)" >> "$P"
  for a in 1 2 3; do
    $BIN/ReLERNN_PREDICT -v "$vcf" -d "$proj" --unphased --batchSizeOverride 50 --gpuID $GPU \
        > "$L/chr42_${sex}_auto.pred.log" 2>&1 && break
    echo "[$sex] PREDICT died attempt $a" >> "$P"; sleep 15
  done
  pred=$proj/chr42_${sex}.PREDICT.txt
  cp -f "$pred" "$OUT/chr42_${sex}.autonet.PREDICT.txt" 2>/dev/null
  if [ -s "$OUT/chr42_${sex}.autonet.PREDICT.txt" ]; then
    awk -v s="$sex" 'NR==1{for(i=1;i<=NF;i++)if($i=="recombRate")r=i;next}
      {rate=$r*1e8;tot+=rate;n++} END{printf "[%s] chr42 via autosomal net: mean=%.4f cM/Mb (n=%d)\n",s,tot/n,n}' \
      "$OUT/chr42_${sex}.autonet.PREDICT.txt" >> "$P"
  else echo "[$sex] PREDICT FAILED" >> "$P"; fi
done
echo "DONE $(date -u +%FT%TZ)" >> "$P"
