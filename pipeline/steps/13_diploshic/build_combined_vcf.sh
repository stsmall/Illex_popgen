#!/usr/bin/env bash
set -euo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh
# baker-633 = callset samples EXCLUDING the 20 small_2026 `ill_*` samples (sex irrelevant for autosomal scan)
$BB/bcftools query -l $ANA/steps/00_callset/filtered/1/variants_filt.vcf.gz | grep -v '^ill_' > $WD/all_samples.txt
[ "$(wc -l < $WD/all_samples.txt)" -eq 633 ] || { echo "ERROR: expected 633 baker samples, got $(wc -l < $WD/all_samples.txt)"; exit 1; }
mkdir -p $WD/combined/parts
do_chr(){ local c=$1
  $BB/bcftools view -S $WD/all_samples.txt --force-samples -m2 -M2 -v snps \
    $ANA/steps/00_callset/filtered/$c/variants_filt.vcf.gz -Ou 2>/dev/null \
  | $BB/bcftools view -e 'F_MISSING>0.5' -Ou 2>/dev/null \
  | $BB/bcftools +fill-tags -Ou -- -t MAF 2>/dev/null \
  | $BB/bcftools view -e 'INFO/MAF<0.01' -Oz -o $WD/combined/all.$c.vcf.gz 2>/dev/null
  $BB/tabix -f -p vcf $WD/combined/all.$c.vcf.gz
}
export -f do_chr; export BB ANA WD
printf '%s\n' $AUT | xargs -P8 -I{} bash -c 'do_chr "$@"' _ {}
echo "done: $(ls $WD/combined/all.*.vcf.gz | wc -l) chrom VCFs"
