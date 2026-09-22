#!/bin/bash
# pi in blocks vs flanks at MAPQ>=20 (as in the paper) and MAPQ>=60 (unique mappers only); 25 individuals (6A subsample).
T=/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_perpop/tests; mkdir -p $T/angsd
ANGSD=/home/ssmall/programs/angsd/angsd; REALSFS=/home/ssmall/programs/angsd/misc/realSFS; TS=/home/ssmall/programs/angsd/misc/thetaStat
ANA=/sietch_colab/data_share/illex/popgen_data/analysis
REF=/sietch_colab/data_share/illex/popgen_data/seq_data/resources/Illex_F24.primary.clean.fa
SITES=$ANA/steps/_shared/accessible_sites.sites
BAMS=$ANA/steps/08_demography/per_population/bamlists_sub/6A.bamlist.txt
grep -E "chr30:|chr23:0.1|chr38:61|chr9:|chr35:|chr1:" $T/intervals.bed > $T/angsd_intervals.bed
one(){ c=$1; s=$2; e=$3; q=$4; o=$T/angsd/${c}_${s}_${e}_q${q}
  $ANGSD -b $BAMS -sites $SITES -r $c:$((s+1))-$e -anc $REF -ref $REF -out $o -doSaf 1 -doCounts 1 -GL 1 -P 2 -minInd 1 -setMinDepthInd 1 -minQ 20 -minMapQ $q -remove_bads 1 -only_proper_pairs 1 >/dev/null 2>&1 || return
  $REALSFS $o.saf.idx -fold 1 -P 2 > $o.sfs 2>/dev/null
  $REALSFS saf2theta $o.saf.idx -sfs $o.sfs -fold 1 -outname $o >/dev/null 2>&1
  $TS do_stat $o.thetas.idx >/dev/null 2>&1; echo "$c $s $e $q $(tail -1 $o.thetas.idx.pestPG | cut -f4,5,9,14)"; }
export -f one; export T ANGSD REALSFS TS REF SITES BAMS
while read c s e blk part; do for q in 20 60; do echo "$c $s $e $q"; done; done < $T/angsd_intervals.bed | xargs -P 8 -n 4 bash -c 'one "$@"' _ > $T/angsd_mq_results.txt
echo ANGSD_DONE
