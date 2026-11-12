#!/bin/bash
# As steps/03_karyotype/inversion_content/extract_af.sh, with the karyotype_final groups (AA223, BB97).
BCF=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
V=/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/2/variants_filt.vcf.gz
L=/sietch_colab/ssmall/projects/msinv_dir/msinv_pboot_wt/results/illex/observed_karyofix/lists
O=/sietch_colab/ssmall/projects/msinv_dir/msinv_pboot_wt/.tmp/karyofix/fst
for g in AA223 BB97; do
  $BCF view -S $L/${g}_samples.txt --force-samples "$V" -Ou 2>/dev/null | $BCF +fill-tags -Ou -- -t AF 2>/dev/null | $BCF query -f '%CHROM\t%POS\t%AF\n' 2>/dev/null > $O/${g}_chr2_af.txt &
done
wait
echo "EXTRACT_DONE $(wc -l < $O/AA223_chr2_af.txt) $(wc -l < $O/BB97_chr2_af.txt)" > $O/extract_status.txt
