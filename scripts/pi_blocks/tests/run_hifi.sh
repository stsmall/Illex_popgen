#!/bin/bash
# HiFi (assembly individual F_24 and male M_71, both aligned to Illex_F24): index, then per-interval depth
# (MAPQ>=0 and >=20), and heterozygous-site density + allele fractions from a pileup over each interval.
T=/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_perpop/tests
H=/sietch_colab/data_share/illex/seq_data/HIFI; S=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/samtools; B=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
REF=/sietch_colab/data_share/illex/popgen_data/seq_data/resources/Illex_F24.primary.clean.fa
F=$(ls $H/F_24/F_24.*.mapped.srt.bam); M=$(ls $H/M_71/*.bam.srt)
echo "index start $(date)"
for f in $F $M; do [ -e $f.bai ] || [ -e $f.csi ] || echo $f; done | xargs -P 4 -I{} $S index -@ 4 -c {}
echo "index done $(date)"
for ind in F_24 M_71; do
  [ $ind = F_24 ] && L="$F" || L="$M"
  for q in 0 20; do $S bedcov -Q $q $T/iv3.bed $L > $T/hifi_${ind}_bedcov_q${q}.tsv; done
  echo "$ind bedcov done $(date)"
  # het sites: mpileup over intervals (skip huge ones by region list), min depth 8, per-allele depth
  $B mpileup -f $REF -R $T/iv3.bed -a AD,DP -q 20 -Q 20 -d 500 --threads 4 -Ou $L 2>/dev/null | \
    $B call -m -v -Ou --threads 4 2>/dev/null | $B view -i 'GT="het" && FMT/DP>=8' -Ou | \
    $B query -f '%CHROM\t%POS\t[%AD]\t[%DP]\n' > $T/hifi_${ind}_hets.tsv
  echo "$ind hets done $(date)"
done
echo HIFI_DONE
