#!/bin/bash
while [ ! -f /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/extract_status.txt ]; do sleep 15; done
cd /sietch_colab/data_share/illex/popgen_data/analysis/manuscript/scripts
env MPLCONFIGDIR=/dev/shm/mplcache /home/ssmall/miniforge3/envs/mkado-vcf/bin/python plot_chr2_karyotype_fst.py > /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/plot_karyotype_fst.log 2>&1
echo "PLOT_DONE $(date)" >> /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/plot_karyotype_fst.log
