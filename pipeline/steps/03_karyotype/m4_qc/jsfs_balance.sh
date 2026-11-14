#!/usr/bin/env bash
# Test whether the private-rare asymmetry of the collinear joint spectrum (jsfs_01 0.352 vs jsfs_10 0.272 for the
# original AA40 x BB40 sets) comes from the sampling-event composition of the two 40-sets (review M4).
# Same ANGSD SAF filters, sites and realSFS 2D call as msinv illex/abc/smc_results/baseline/run_jsfs_control2.sh,
# on one 10-Mb collinear block (chr2:85-95 Mb), for the original sets and three division-balanced draws
# (identical division counts in both groups, lists drawn with seeds 20261114-20261116, in
# analysis/steps/03_karyotype/m4_qc/jsfs_balance). realSFS capped at 100 iterations. Cells: jsfs_balance.py.
set -uo pipefail
A=/sietch_colab/data_share/illex/popgen_data/analysis
O=$A/steps/03_karyotype/m4_qc/jsfs_balance; cd "$O"
REF=/sietch_colab/data_share/illex/andy_links/Illex_F24.primary.clean.fa
ANGSD=/home/ssmall/programs/angsd/angsd
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
SITES=$A/steps/04_angsd_chr2/chr2.accessible.sites
REG=2:85000001-95000000
saf(){ local g=$1 s=$2; [ -s "$g.$s.saf.idx" ] && return
  "$ANGSD" -b "${g}40_$s.txt" -sites "$SITES" -r "$REG" -anc "$REF" -ref "$REF" -out "$g.$s" \
    -doSaf 1 -doCounts 1 -GL 1 -P 4 -minInd 1 -setMinDepthInd 1 -minQ 20 -minMapQ 20 \
    -remove_bads 1 -only_proper_pairs 1 > "$g.$s.angsd.log" 2>&1; }
export -f saf; export ANGSD REF SITES REG
for s in orig bal1 bal2 bal3; do for g in AA BB; do echo "$g $s"; done; done | xargs -P 4 -n 2 bash -c 'saf "$@"' _
echo "[$(date +%T)] SAFs done"
sfs2d(){ local s=$1; [ -s "AA_BB.$s.2dSFS" ] && return
  "$REALSFS" "AA.$s.saf.idx" "BB.$s.saf.idx" -fold 0 -P 4 -tole 1e-8 -maxIter 100 \
    > "AA_BB.$s.2dSFS.part" 2> "AA_BB.$s.log" && mv "AA_BB.$s.2dSFS.part" "AA_BB.$s.2dSFS"; }
export -f sfs2d; export REALSFS
printf '%s\n' orig bal1 bal2 bal3 | xargs -P 4 -I{} bash -c 'sfs2d "$@"' _ {}
echo "[$(date +%T)] 2dSFS done"
touch ALL_DONE
