#!/bin/bash
cd /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/flank_div
BCF=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
V=/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/2/variants_filt.vcf.gz
P=/sietch_colab/ssmall/projects/msinv_dir/inversion_sims/files/.tmp/illex_chr2/pops.tsv
cut -f1 $P > karyo_samples.txt
echo "extract start $(date)"
$BCF view -r 2:40000000-100000000 -S karyo_samples.txt --force-samples --threads 8 -Oz -o flank_40_100.vcf.gz "$V" && $BCF index -t flank_40_100.vcf.gz
echo "extract done $(date) $($BCF index -n flank_40_100.vcf.gz) variants"
CUDA_VISIBLE_DEVICES=0 /home/ssmall/miniforge3/envs/varbuddy-pggpu/bin/python compute_flank_div2.py
echo "V2_DONE $(date)"
