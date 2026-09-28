#!/usr/bin/env bash
# chr2 3-state accessibility+called mask for msinv:
#   accessible_variant | accessible_invariant | inaccessible  (tiles chr2 completely)
# From accessible_sites.bed (callable footprint) + all.2.vcf.gz (variant positions).
set -uo pipefail
B=/sietch_colab/data_share/illex/popgen_data
FAI=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa.fai
ACC=$B/degenotate_illex/accessible_sites.bed
VCF=$B/seq_data/baker_2025/vcfs/all.2.vcf.gz
BT=/home/ssmall/miniforge3/envs/annotation-expression-buddy-barrnap_trnascan/bin/bedtools
BCF=/home/ssmall/bin/bcftools
OUT=$B/analysis/steps/03_karyotype/chr2_mask; mkdir -p "$OUT"; cd "$OUT"
export LC_ALL=C
echo "MASK_START $(date -u +%FT%TZ)" > progress
L=$(awk '$1=="2"{print $2}' "$FAI"); printf '2\t%s\n' "$L" > genome2.txt
# accessible chr2 (sorted + merged)
awk '$1=="2"' "$ACC" | sort -k1,1 -k2,2n | "$BT" merge -i - > acc2.bed
# chr2 variant positions -> 0-based BED (threaded decompress of the 42GB VCF)
"$BCF" query -f '%CHROM\t%POS\n' "$VCF" 2>>progress \
  | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' | sort -k1,1 -k2,2n > var2.bed
nvar=$(wc -l < var2.bed); echo "VAR_DONE $nvar $(date -u +%FT%TZ)" >> progress
# 3-state partition (non-overlapping, complete tiling)
"$BT" intersect -a var2.bed -b acc2.bed -sorted | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"accessible_variant"}'   > s_av.bed
"$BT" subtract  -a acc2.bed -b var2.bed              | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"accessible_invariant"}' > s_ai.bed
"$BT" complement -i acc2.bed -g genome2.txt          | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"inaccessible"}'         > s_in.bed
cat s_av.bed s_ai.bed s_in.bed | sort -k1,1 -k2,2n > chr2.mask.3state.bed
{ echo "=== chr2.mask.3state.bed summary ==="
  awk '{bp=$3-$2; s[$4]+=bp; n[$4]++} END{for(k in s) printf "  %-22s %11d bp  %9d intervals\n",k,s[k],n[k]}' chr2.mask.3state.bed
  awk -v L="$L" '{t+=$3-$2} END{printf "  total tiled: %d bp (chr2 len %d; match=%s)\n",t,L,(t==L?"YES":"NO")}' chr2.mask.3state.bed
  echo "  n variants: $nvar ; accessible_variant sites: $(wc -l < s_av.bed)"
} | tee -a progress
echo "MASK_DONE $(date -u +%FT%TZ)" >> progress
