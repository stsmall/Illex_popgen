#!/usr/bin/env bash
# Codon-aware Illex vs Sthenoteuthis alignment per gene with MACSE (Sthenoteuthis = less reliable,
# frameshifts allowed). 32 genes at a time (shared machine, <=100 cores in total).
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/18_mutation_bias/stheno/pairs
M=/home/ssmall/miniforge3/envs/mkado_env/bin/macse
export D M
awk 'NR>1{print $1}' $D/pairs.tsv | xargs -P 32 -I{} bash -c '[ -s $D/{}.NT.fa ] || timeout 1800 $M -prog alignSequences -seq $D/{}.ill.fa -seq_lr $D/{}.sth.fa -out_NT $D/{}.NT.fa -out_AA $D/{}.AA.fa > /dev/null 2>&1 || echo FAIL {}'
touch $D/MACSE_DONE
