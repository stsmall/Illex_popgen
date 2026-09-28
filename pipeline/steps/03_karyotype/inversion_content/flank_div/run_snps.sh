#!/bin/bash
cd /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/flank_div
BCF=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
echo "filter start $(date)"
$BCF view -m2 -M2 -v snps --threads 8 flank_40_100.vcf.gz -Ou | $BCF annotate -x 'INFO,^FORMAT/GT' --threads 4 -Oz -o flank_40_100.snps.vcf.gz && $BCF index -t flank_40_100.snps.vcf.gz
echo "filter done $(date) $($BCF index -n flank_40_100.snps.vcf.gz) biallelic SNPs"
CUDA_VISIBLE_DEVICES=0 /home/ssmall/miniforge3/envs/varbuddy-pggpu/bin/python compute_flank_div_snps.py
echo "SNPS_DONE $(date)"
