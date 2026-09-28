#!/usr/bin/env bash
# Which arrangement (AA vs BB) is ANCESTRAL? At chr2:60-80Mb arrangement-diagnostic
# SNPs (near-fixed differences between AA-homs and BB-homs), tally which arrangement
# carries the est-sfs ancestral allele. Cross-check with coindetii outgroup allele.
set -euo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
BCF=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
K=$ANA/steps/03_karyotype
VCF=$ANA/steps/00_callset/filtered/2/variants_filt.vcf.gz
POL=/sietch_colab/data_share/illex/popgen_data/polarize/illecebrosus_input/gatk_vcfs/illecebrosus.polarized.vcf.gz
COIN=/sietch_colab/data_share/illex/popgen_data/mkado_illex/inputs/coindetii.vcf.gz
REG=2:60000000-80000000
OUT=$K/polarize; mkdir -p "$OUT"

echo "[$(date)] AF per arrangement in $REG"
$BCF view -r $REG -S "$K/AA_samples.txt" --force-samples "$VCF" -Ou \
  | $BCF +fill-tags -Ou -- -t AF | $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > "$OUT/aa_af.txt"
$BCF view -r $REG -S "$K/BB_samples.txt" --force-samples "$VCF" -Ou \
  | $BCF +fill-tags -Ou -- -t AF | $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > "$OUT/bb_af.txt"
echo "[$(date)] ancestral (est-sfs AA=) + coindetii allele"
$BCF query -r $REG -f '%CHROM\t%POS\t%INFO/AA\n' "$POL" 2>/dev/null > "$OUT/anc_estsfs.txt"
$BCF query -r $REG -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$COIN" 2>/dev/null > "$OUT/coindetii.txt"
echo "  aa_af=$(wc -l <"$OUT/aa_af.txt")  bb_af=$(wc -l <"$OUT/bb_af.txt")  est-sfs anc=$(wc -l <"$OUT/anc_estsfs.txt")  coindetii=$(wc -l <"$OUT/coindetii.txt")"
echo "[$(date)] done extraction -> run polarize_arrangements.py"
