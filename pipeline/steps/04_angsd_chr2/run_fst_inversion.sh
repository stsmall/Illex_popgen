#!/usr/bin/env bash
# FST between karyotype arrangements, RESTRICTED to the inversion region 2:60-80Mb.
# Reuses existing whole-chr2 SAFs; realSFS -r keeps the 2D-SFS EM small & stable.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
D=$ANA/steps/04_angsd_chr2; SAF=$D/saf; F=$D/fst; mkdir -p "$F"
REG=2:60000000-80000000; TH=8
for pair in "AA BB" "AA AB" "AB BB"; do
  set -- $pair; p1=$1; p2=$2; o="$F/${p1}_${p2}.inv"
  echo "[$(date)] 2dSFS $p1 vs $p2 ($REG)"
  "$REALSFS" "$SAF/${p1}.2.saf.idx" "$SAF/${p2}.2.saf.idx" -r "$REG" -fold 1 -P $TH -tole 1e-8 -maxIter 200 > "$o.2dSFS" 2>"$D/logs/fstinv.${p1}_${p2}.log"
  "$REALSFS" fst index "$SAF/${p1}.2.saf.idx" "$SAF/${p2}.2.saf.idx" -r "$REG" -sfs "$o.2dSFS" -fstout "$o" -fold 1 2>>"$D/logs/fstinv.${p1}_${p2}.log"
  "$REALSFS" fst stats "$o.fst.idx" > "$o.fst.txt" 2>>"$D/logs/fstinv.${p1}_${p2}.log"
  echo "  FST $p1 vs $p2: $(cat "$o.fst.txt")"
done
echo "[$(date)] FST(inversion) done"
