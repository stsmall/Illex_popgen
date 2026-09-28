#!/bin/bash
echo -e "sample\tchr1_block\tchr1_ctrl\tchr34_block\tchr34_ctrl" > /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_depth.tsv
ls /sietch_colab/data_share/illex/popgen_data/seq_data/baker_2025/mapping/dedup/*.bam | xargs -P 16 -n 1 bash /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_depth_one.sh >> /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_depth.tsv
echo DEPTH_DONE >> /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_depth.tsv
