#!/bin/bash
cd /sietch_colab/data_share/illex/popgen_data/analysis/steps/_shared/genmap
B=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bedtools
ACC=/sietch_colab/data_share/illex/popgen_data/degenotate_illex/accessible_sites.bed
BLK=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/config/multicopy_block_mask.bed
awk 'BEGIN{OFS="\t"} $4<1 {print $1,$2,$3}' illex_k150_e2.bedgraph | $B merge -i - > nonunique_k150_e2.bed
awk 'BEGIN{OFS="\t"} $4<=0.5 {print $1,$2,$3}' illex_k150_e2.bedgraph | $B merge -i - > low_k150_e2.bed
tot=$(awk '{s+=$3-$2} END{print s}' illex_k150_e2.bedgraph)
nu=$(awk '{s+=$3-$2} END{print s}' nonunique_k150_e2.bed); lo=$(awk '{s+=$3-$2} END{print s}' low_k150_e2.bed)
acc=$(awk '{s+=$3-$2} END{print s}' $ACC)
sort -k1,1 -k2,2n $ACC > acc.sorted.bed; sort -k1,1 -k2,2n $BLK > blk.sorted.bed
nuacc=$($B intersect -a nonunique_k150_e2.bed -b acc.sorted.bed | awk '{s+=$3-$2} END{print s+0}')
loacc=$($B intersect -a low_k150_e2.bed -b acc.sorted.bed | awk '{s+=$3-$2} END{print s+0}')
blk=$(awk '{s+=$3-$2} END{print s}' blk.sorted.bed)
nublk=$($B intersect -a nonunique_k150_e2.bed -b blk.sorted.bed | awk '{s+=$3-$2} END{print s+0}')
accblk=$($B intersect -a acc.sorted.bed -b blk.sorted.bed | awk '{s+=$3-$2} END{print s+0}')
nuaccblk=$($B intersect -a nonunique_k150_e2.bed -b acc.sorted.bed | $B intersect -a - -b blk.sorted.bed | awk '{s+=$3-$2} END{print s+0}')
echo -e "genome_bp\t$tot\nnonunique_bp\t$nu\nlow_le0.5_bp\t$lo\naccessible_bp\t$acc\nnonunique_in_accessible_bp\t$nuacc\nlow_in_accessible_bp\t$loacc\nblocks_bp\t$blk\nnonunique_in_blocks_bp\t$nublk\naccessible_in_blocks_bp\t$accblk\nnonunique_accessible_in_blocks_bp\t$nuaccblk" > genmap_summary.tsv
echo SUMMARY_DONE
