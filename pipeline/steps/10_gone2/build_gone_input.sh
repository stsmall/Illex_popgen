#!/usr/bin/env bash
# Build GONE2 input: baker-633, biallelic SNPs, F_MISSING<0.2, MAF>=0.05, autosomes (1,3-41,43-45;
# excl chr2 inversion + chr42 sex + chrZ), accessible sites (variants_filt already accessible-masked).
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
GD=$A/steps/10_gone2; mkdir -p $GD/parts
AUT="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"
export BB A GD
do_chr(){ local c=$1
  $BB/bcftools view -S $A/steps/00_callset/baker_633.txt --force-samples -m2 -M2 -v snps \
    $A/steps/00_callset/filtered/$c/variants_filt.vcf.gz -Ou 2>/dev/null \
  | $BB/bcftools view -e 'F_MISSING>0.2' -Ou 2>/dev/null \
  | $BB/bcftools +fill-tags -Ou -- -t MAF 2>/dev/null \
  | $BB/bcftools view -e 'INFO/MAF<0.05' -Oz -o $GD/parts/$c.vcf.gz 2>/dev/null
  $BB/tabix -f -p vcf $GD/parts/$c.vcf.gz 2>/dev/null
  echo "chr$c: $($BB/bcftools index -n $GD/parts/$c.vcf.gz) SNPs"
}
export -f do_chr
echo "[$(date)] building GONE2 parts (-P8)"
printf '%s\n' $AUT | xargs -P8 -I{} bash -c 'do_chr "$@"' _ {}
echo "[$(date)] concat"
ls $GD/parts/*.vcf.gz | sort -t/ -k999 > /dev/null
$BB/bcftools concat -Oz -o $GD/gone_input.vcf.gz $(for c in $AUT; do echo $GD/parts/$c.vcf.gz; done) 2>/dev/null
$BB/tabix -f -p vcf $GD/gone_input.vcf.gz 2>/dev/null
echo "[$(date)] TOTAL SNPs: $($BB/bcftools index -n $GD/gone_input.vcf.gz)  samples: $($BB/bcftools query -l $GD/gone_input.vcf.gz|wc -l)"
