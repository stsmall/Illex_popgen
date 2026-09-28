#!/usr/bin/env bash
# Dense per-sex autosomal VCFs for ReLERNN: biallelic SNPs, F_MISSING<0.5, MAF>0.01 (denser than
# the 0.05 GONE set -> finer windows), autosomes (excl chr2 inv, chr42, chrZ). Missing kept as './.'.
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
SD=$A/steps/12_sex; RD=$A/steps/11_relernn
AUT="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"
sex=$1; slist=$SD/${sex}s.txt          # males.txt / females.txt (confident)
outdir=$RD/$sex; mkdir -p $outdir/parts
export BB A slist outdir
do_chr(){ local c=$1
  $BB/bcftools view -S "$slist" --force-samples -m2 -M2 -v snps \
    $A/steps/00_callset/filtered/$c/variants_filt.vcf.gz -Ou 2>/dev/null \
  | $BB/bcftools view -e 'F_MISSING>0.5' -Ou 2>/dev/null \
  | $BB/bcftools +fill-tags -Ou -- -t MAF 2>/dev/null \
  | $BB/bcftools view -e 'INFO/MAF<0.01' -Oz -o $outdir/parts/$c.vcf.gz 2>/dev/null
  $BB/tabix -f -p vcf $outdir/parts/$c.vcf.gz 2>/dev/null
}
export -f do_chr
printf '%s\n' $AUT | xargs -P8 -I{} bash -c 'do_chr "$@"' _ {}
$BB/bcftools concat -Oz -o $outdir/${sex}.vcf.gz $(for c in $AUT; do echo $outdir/parts/$c.vcf.gz; done) 2>/dev/null
$BB/bgzip -dk -f $outdir/${sex}.vcf.gz
echo "[$(date)] $sex: $($BB/bcftools index -n $outdir/${sex}.vcf.gz 2>/dev/null || grep -vc '^#' $outdir/${sex}.vcf) SNPs, $($BB/bcftools query -l $outdir/${sex}.vcf.gz|wc -l) samples"
