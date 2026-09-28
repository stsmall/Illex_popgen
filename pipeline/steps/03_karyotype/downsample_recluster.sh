#!/usr/bin/env bash
# Depth-control for chr2 inversion karyotyping: downsample all baker BAMs to a
# common ~1x depth in the inversion region, re-genotype, re-cluster. If the 3
# karyotype clusters persist at UNIFORM depth, the inversion is not a depth artifact.
# SHARED MACHINE: cap parallelism (<=16 downsample jobs; mpileup 8 threads) -> <=100 cores.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy
SAM="$BB/bin/samtools"; BCF="$BB/bin/bcftools"
REGION="2:60000000-80000000"; REGLEN=20000000; TARGET=1.0; SEED=42
DEDUP=/sietch_colab/data_share/illex/popgen_data/seq_data/baker_2025/mapping/dedup
DS=$ANA/steps/03_karyotype/downsample
mkdir -p "$DS/bams" "$DS/logs"
BAKER="$ANA/steps/00_callset/baker_633.txt"

echo "[$(date)] STEP1+2: region depth + downsample to ${TARGET}x"
downsample_one() {
  local s=$1
  local bam="$DEDUP/$s.bam"
  local out="$DS/bams/$s.ds.bam"
  [ -f "$out.bai" ] && return 0
  [ -f "$bam" ] || { echo "MISSING $bam"; return 0; }
  local reads=$("$SAM" view -c "$bam" "$REGION" 2>/dev/null)
  local dep=$(awk -v r="$reads" -v L="$REGLEN" 'BEGIN{print r*150/L}')
  local frac=$(awk -v t="$TARGET" -v d="$dep" 'BEGIN{f=(d>0)?t/d:1; if(f>=1)f=1; printf "%.4f",f}')
  if awk -v f="$frac" 'BEGIN{exit !(f>=1)}'; then
    "$SAM" view -b "$bam" "$REGION" > "$out" 2>/dev/null   # already <=target: keep all
  else
    "$SAM" view -s "${SEED}${frac#0}" -b "$bam" "$REGION" > "$out" 2>/dev/null  # -s SEED.FRAC
  fi
  "$SAM" index "$out"
}
export -f downsample_one; export SAM DS DEDUP REGION REGLEN TARGET SEED
tail -n +1 "$BAKER" | xargs -P 16 -I{} bash -c 'downsample_one "$@"' _ {} > "$DS/logs/downsample.log" 2>&1

ls "$DS/bams/"*.ds.bam > "$DS/ds_bamlist.txt"
echo "[$(date)] downsampled bams: $(wc -l < "$DS/ds_bamlist.txt")"

echo "[$(date)] STEP3: mpileup+call on depth-equalized BAMs (region ${REGION})"
DSVCF="$DS/chr2inv_ds.vcf.gz"
"$BCF" mpileup -b "$DS/ds_bamlist.txt" -r "$REGION" -f "$REF" -a FORMAT/DP --threads 8 -q 20 -Q 20 2>"$DS/logs/mpileup.log" \
  | "$BCF" call -mv --threads 8 2>>"$DS/logs/mpileup.log" \
  | "$BCF" view -m2 -M2 -v snps -q 0.05:minor -i 'F_MISSING<0.5' -Oz -o "$DSVCF" 2>>"$DS/logs/mpileup.log"
"$BCF" index -t "$DSVCF"
echo "[$(date)] downsampled SNP VCF: $("$BCF" index -n "$DSVCF") biallelic MAF>0.05 SNPs"
echo "[$(date)] DONE downsample+call. Next: pg_gpu PCA + GMM (recluster_downsampled.py)"
