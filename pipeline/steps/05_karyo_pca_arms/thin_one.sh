#!/usr/bin/env bash
set -uo pipefail
c=$1
ANA=/sietch_colab/data_share/illex/popgen_data/analysis
BBIN=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
BAKER=$ANA/steps/00_callset/baker_633.txt
TH=$ANA/steps/05_karyo_pca_arms/thinned
out=$TH/chr${c}.thin.vcf.gz
[[ -s "$out.tbi" ]] && { echo "chr$c exists"; exit 0; }
"$BBIN/bcftools" view -S "$BAKER" --force-samples -m2 -M2 -v snps -q 0.05:minor \
    "$ANA/steps/00_callset/filtered/${c}/variants_filt.vcf.gz" -Ou 2>/dev/null \
  | "$BBIN/bcftools" +prune -w 10kb -n 1 -Oz -o "$out" 2>/dev/null
"$BBIN/bcftools" index -t "$out" && echo "chr$c: $("$BBIN/bcftools" index -n "$out") SNPs" || echo "chr$c FAILED"
