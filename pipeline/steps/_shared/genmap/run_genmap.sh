#!/bin/bash
# GenMap mappability of the Illex_F24 assembly for 150-bp reads (k=150, up to 2 mismatches).
G=/home/ssmall/miniforge3/envs/assembly-buddy-mappability/bin/genmap
REF=/sietch_colab/data_share/illex/popgen_data/seq_data/resources/Illex_F24.primary.clean.fa
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/_shared/genmap
echo "index start $(date)"
[ -d $D/index ] || $G index -F $REF -I $D/index
echo "index done $(date)"
$G map -K 150 -E 2 -I $D/index -O $D/illex_k150_e2 -bg -T 32
echo "GENMAP_DONE $(date)"
