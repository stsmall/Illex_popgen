#!/usr/bin/env bash
# chrZ diversity by sex via ANGSD SAF->realSFS->theta. Males (ZZ, diploid) vs Females (ZW).
# Accessible chrZ sites (genome-wide .bin), folded. Compare pi_Z to autosomal pi (0.0093).
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
source $A/config.sh
ANGSD=/home/ssmall/programs/angsd/angsd
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
THETASTAT=/home/ssmall/programs/angsd/misc/thetaStat
SD=$A/steps/12_sex; SAF=$SD/z_saf; TH=$SD/z_thetas; mkdir -p $SAF $TH $SD/logs
SITES=$A/steps/_shared/accessible_sites.sites   # genome-wide indexed (.bin/.idx)
touch "$REF.fai"
do_sex(){ local sex=$1 bl=$2
  $ANGSD -b "$bl" -sites "$SITES" -r Z -anc "$REF" -ref "$REF" -out "$SAF/z_$sex" \
    -doSaf 1 -doCounts 1 -GL 1 -P 6 -minInd 1 -setMinDepthInd 1 \
    -minQ 20 -minMapQ 20 -remove_bads 1 -only_proper_pairs 1 > "$SD/logs/z_$sex.log" 2>&1
  $REALSFS "$SAF/z_$sex.saf.idx" -fold 1 -P 6 -tole 1e-8 -maxIter 200 > "$TH/z_$sex.sfs" 2>>"$SD/logs/z_$sex.log"
  $REALSFS saf2theta "$SAF/z_$sex.saf.idx" -sfs "$TH/z_$sex.sfs" -outname "$TH/z_$sex" -fold 1 2>>"$SD/logs/z_$sex.log"
  $THETASTAT do_stat "$TH/z_$sex.thetas.idx" -outnames "$TH/z_$sex.avg" 2>>"$SD/logs/z_$sex.log"
  echo "[$(date)] chrZ $sex done: $(cat $TH/z_$sex.avg.pestPG 2>/dev/null | sed -n 2p | awk '{printf "tW=%.0f tP=%.0f TajD=%.3f nSites=%d pi=%.6f", $4,$5,$9,$14,$5/$14}')"
}
echo "[$(date)] chrZ by sex: males (n=$(wc -l<$SD/males_bamlist.txt)) + females (n=$(wc -l<$SD/females_bamlist.txt))"
do_sex male   $SD/males_bamlist.txt &
do_sex female $SD/females_bamlist.txt &
wait
echo "[$(date)] chrZ by-sex DONE"
