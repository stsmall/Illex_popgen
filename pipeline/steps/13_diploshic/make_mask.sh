#!/usr/bin/env bash
set -euo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh
BED=/sietch_colab/data_share/illex/popgen_data/degenotate_illex/accessible_sites.bed  # verified path

# BUG FIX: the accessible BED chroms are LEXICALLY ordered (1,10,11,..,2,20,..) but ${REF}.fai is
# NUMERIC (1,2,3,..). `bedtools complement` requires -i in the SAME chrom order as -g, else it emits
# whole chromosomes (3-9) as complement -> mask N's them 100% -> fvecVcf gets 0 SNPs. Sort BED to fai order.
$BB/bedtools sort -i $BED -g ${REF}.fai > $WD/accessible.faidxsorted.bed
$BB/bedtools complement -i $WD/accessible.faidxsorted.bed -g ${REF}.fai > $WD/inaccessible.bed
$BB/bedtools maskfasta -fi $REF -bed $WD/inaccessible.bed -fo $WD/mask.fa
$BB/samtools faidx $WD/mask.fa
# SANITY: chr3-9 (the ones the bug wiped) must NOT be ~100% N
for c in 3 5 9; do
  pctN=$($BB/samtools faidx $WD/mask.fa $c 2>/dev/null | tail -n +2 | tr -d '\n' | awk '{t=length($0);g=gsub(/N/,"");printf "%.1f", (t? g*100.0/t : 100)}')
  echo "chr$c: ${pctN}% N"
  awk -v p="$pctN" 'BEGIN{if(p>99){print "ERROR: chr sample still ~100% N — mask fix FAILED"; exit 1}}'
done
echo "mask.fa built and indexed (chr3-9 accessibility restored)"
