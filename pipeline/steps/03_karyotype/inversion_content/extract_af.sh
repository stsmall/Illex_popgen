#!/bin/bash
BCF=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
V=/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/2/variants_filt.vcf.gz
K=/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype
O=/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content
$BCF view -S $K/AA_samples.txt --force-samples "$V" -Ou 2>/dev/null | $BCF +fill-tags -Ou -- -t AF 2>/dev/null | $BCF query -f '%CHROM\t%POS\t%AF\n' 2>/dev/null > $O/aa_chr2_af.txt
$BCF view -S $K/BB_samples.txt --force-samples "$V" -Ou 2>/dev/null | $BCF +fill-tags -Ou -- -t AF 2>/dev/null | $BCF query -f '%CHROM\t%POS\t%AF\n' 2>/dev/null > $O/bb_chr2_af.txt
echo "EXTRACT_DONE $(wc -l < $O/aa_chr2_af.txt) $(wc -l < $O/bb_chr2_af.txt)" > $O/extract_status.txt
