#!/bin/bash
# HiFi tests, part 2: M_71 per-interval depth, and heterozygous-site density/allele fractions in 1-Mb samples
# of key blocks and flanks for both individuals (F_24 = the assembly individual; M_71 = a second individual).
T=/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_perpop/tests
H=/sietch_colab/data_share/illex/seq_data/HIFI
S=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/samtools
B=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
REF=/sietch_colab/data_share/illex/popgen_data/seq_data/resources/Illex_F24.primary.clean.fa
M=$(echo $(ls $H/M_71/*.bam.srt))
F=$(echo $(ls $H/F_24/F_24.*.mapped.srt.bam))
if ! grep -q "M_71 bedcov done" $T/run_hifi2.log 2>/dev/null && ! pgrep -f "bedcov -Q 0 $T/iv3.bed $H/M_71" >/dev/null; then
( for q in 0 20; do $S bedcov -Q $q $T/iv3.bed $M > $T/hifi_M_71_bedcov_q${q}.tsv; done; echo "M_71 bedcov done $(date)" ) &
fi
het() {
  name=$1; c=$2; s=$3; e=$4; ind=$5; shift 5
  tmp=$T/tmp_${ind}_${name}.bam; $S merge -f -R $c:$s-$e -o $tmp "$@" 2>/dev/null && $S index $tmp
  $B mpileup -f $REF -r $c:$s-$e -a AD,DP -q 20 -Q 20 -d 300 --ignore-RG -Ou $tmp 2>/dev/null | $B call -m -v -Ou 2>/dev/null | \
    $B view -i 'GT="het" && FMT/DP>=8' -Ou | $B query -f "$name\t$ind\t%POS\t[%AD]\t[%DP]\n" > $T/het_${ind}_${name}.tsv
  rm -f $tmp $tmp.bai; echo "het $ind $name done $(date)"
}
N=0
while read n c s e; do
  het $n $c $s $e F_24 $F &
  het $n $c $s $e M_71 $M &
  N=$((N+2)); if [ $N -ge 10 ]; then wait; N=0; fi
done < $T/het_regions.txt
wait
echo HETS_DONE
