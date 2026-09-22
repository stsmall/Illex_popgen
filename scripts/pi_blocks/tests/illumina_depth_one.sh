#!/bin/bash
b=$1; s=$(basename $b .bam); S=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/samtools; T=/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_perpop/tests
q0=$($S bedcov -Q 0 $T/iv3.bed $b 2>/dev/null | cut -f4 | paste -sd,); q20=$($S bedcov -Q 20 $T/iv3.bed $b 2>/dev/null | cut -f4 | paste -sd,)
echo -e "$s\t$q0\t$q20"
