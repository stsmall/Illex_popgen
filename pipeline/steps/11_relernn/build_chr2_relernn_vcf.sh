#!/usr/bin/env bash
# chr2 ReLERNN input VCFs (same filters as build_persex_vcf.sh: biallelic SNP, F_MISSING<0.5, MAF>0.01):
#   AA  = homokaryotype A (n=254) -> within-arrangement (collinear) background rho for msinv
#   all = all samples (mixed karyotype) -> observed inversion-suppressed landscape
set -uo pipefail
B=/sietch_colab/data_share/illex/popgen_data
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
SRC=$B/analysis/steps/00_callset/filtered/2/variants_filt.vcf.gz
OUT=$B/analysis/steps/11_relernn/chr2; mkdir -p "$OUT"
AA=$B/analysis/steps/03_karyotype/AA_samples.txt
build(){ local tag=$1 slist=$2 sopt=""
  [ "$slist" != "ALL" ] && sopt="-S $slist --force-samples"
  $BB/bcftools view $sopt -m2 -M2 -v snps "$SRC" -Ou 2>/dev/null \
    | $BB/bcftools view -e 'F_MISSING>0.5' -Ou 2>/dev/null \
    | $BB/bcftools +fill-tags -Ou -- -t MAF 2>/dev/null \
    | $BB/bcftools view -e 'INFO/MAF<0.01' -Oz -o "$OUT/$tag.chr2.vcf.gz" 2>/dev/null
  $BB/tabix -f -p vcf "$OUT/$tag.chr2.vcf.gz" 2>/dev/null
  $BB/bgzip -dk -f "$OUT/$tag.chr2.vcf.gz"
  echo "[$tag] $($BB/bcftools index -n "$OUT/$tag.chr2.vcf.gz" 2>/dev/null) SNPs, $($BB/bcftools query -l "$OUT/$tag.chr2.vcf.gz" | wc -l) samples" >> "$OUT/build.log"
}
echo "VCFBUILD_START $(date -u +%FT%TZ)" > "$OUT/build.log"
build AA "$AA"
build all ALL
echo "VCFBUILD_DONE $(date -u +%FT%TZ)" >> "$OUT/build.log"
