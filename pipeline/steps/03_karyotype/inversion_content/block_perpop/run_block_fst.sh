#!/bin/bash
cd /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_perpop
BCF=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
V=/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered
TAGS=$(cut -f2 sample_groups.txt | sort -u | awk '{printf "%%AC_%s\\t%%AN_%s\\t",$1,$1}')
for spec in "1:19000000-36300000" "30:12300000-23500000" "23:1-10600000" "38:56500000-73000000"; do
  c=${spec%%:*}; echo "start $spec $(date)"
  $BCF view -r "$spec" -m2 -M2 -v snps --threads 4 "$V/$c/variants_filt.vcf.gz" -Ou | $BCF +fill-tags -Ou -- -S sample_groups.txt -t AC,AN | $BCF query -f "%POS\t${TAGS}\n" > af_chr$c.tsv
  echo "done $spec $(wc -l < af_chr$c.tsv) sites $(date)"
done
echo ALL_DONE
