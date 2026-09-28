#!/usr/bin/env bash
# Per-chrom thinned VCFs for the karyotype-colored PCA panel (chr2-specificity test).
# baker-633 only, biallelic SNPs, MAF>0.05, LD-thin ~1 SNP/10kb (bcftools +prune).
# SHARED MACHINE: -P10.
set -uo pipefail
ANA=/sietch_colab/data_share/illex/popgen_data/analysis
BBIN=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
BAKER=$ANA/steps/00_callset/baker_633.txt
D=$ANA/steps/05_karyo_pca_arms; TH=$D/thinned; mkdir -p "$TH" "$D/logs"
CHRS="2 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 Z"
thin_one() {
  local c=$1 out="$TH/chr${c}.thin.vcf.gz"
  [[ -s "$out.tbi" ]] && return 0
  "$BBIN/bcftools" view -S "$BAKER" --force-samples -m2 -M2 -v snps -q 0.05:minor \
      "$ANA/steps/00_callset/filtered/${c}/variants_filt.vcf.gz" -Ou 2>/dev/null \
    | "$BBIN/bcftools" +prune -w 10000 -n 1 -Oz -o "$out" 2>/dev/null
  "$BBIN/bcftools" index -t "$out"
  echo "chr$c: $("$BBIN/bcftools" index -n "$out") thinned SNPs"
}
export -f thin_one; export BBIN BAKER TH ANA
printf '%s\n' $CHRS | xargs -P 10 -I{} bash -c 'thin_one "$@"' _ {}
echo "[$(date)] thinning done"
