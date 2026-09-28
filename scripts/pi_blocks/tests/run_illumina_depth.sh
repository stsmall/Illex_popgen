#!/bin/bash
T=/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_perpop/tests
ids=/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/baker_633.txt
D=/sietch_colab/data_share/illex/popgen_data/seq_data/baker_2025/mapping/dedup
sed "s#^#$D/#; s#\$#.bam#" $ids | xargs -P 16 -n 1 $T/illumina_depth_one.sh > $T/illumina_depth.tsv
echo ILLUMINA_DONE >> $T/illumina_depth.log
