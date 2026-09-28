#!/usr/bin/env bash
# Idempotent finalize/cleanup for genome-wide baker-633 SAF->realSFS->theta.
# Uses the CANONICAL genome-wide ANGSD sites index (_shared/accessible_sites.sites
# + .bin/.idx built from accessible_sites.bed) via -sites + -r <chrom>.
# Skips chroms whose theta (avg.pestPG) is already done; redoes broken/missing.
# SHARED MACHINE: 4 chroms x -P6 ~= 24 threads. Run AFTER the first sweep frees up.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
ANGSD=/home/ssmall/programs/angsd/angsd
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
THETASTAT=/home/ssmall/programs/angsd/misc/thetaStat
D=$ANA/steps/08_demography
BAM=$ANA/steps/04_angsd_chr2/bamlists/bakerMeta_bamlist.txt   # 633 baker dedup
SITES=$ANA/steps/_shared/accessible_sites.sites               # genome-wide, indexed (.bin/.idx)
SAF=$D/saf; TH=$D/thetas; mkdir -p "$SAF" "$TH" "$D/logs"
touch "$REF.fai"
THREADS=6

do_chrom() {
  local c=$1                       # MUST be its own line: `local c=$1 out=...${c}` expands ${c} BEFORE c is set
  local out="$SAF/baker.${c}"
  [[ -s "$TH/baker.${c}.avg.pestPG" ]] && { echo "[$c] theta done, skip"; return 0; }
  echo "[$(date)] [$c] SAF (633 baker, genome-wide .bin sites, folded)"
  "$ANGSD" -b "$BAM" -sites "$SITES" -r "$c" -anc "$REF" -ref "$REF" -out "$out" \
    -doSaf 1 -doCounts 1 -GL 1 -P "$THREADS" -minInd 1 -setMinDepthInd 1 \
    -minQ 20 -minMapQ 20 -remove_bads 1 -only_proper_pairs 1 > "$D/logs/saf.${c}.log" 2>&1
  [[ -s "$out.saf.idx" ]] || { echo "[$(date)] [$c] FAIL: empty saf.idx"; return 1; }
  "$REALSFS" "$out.saf.idx" -fold 1 -P "$THREADS" -tole 1e-8 -maxIter 200 > "$TH/baker.${c}.sfs" 2>>"$D/logs/saf.${c}.log"
  [[ -s "$TH/baker.${c}.sfs" ]] || { echo "[$(date)] [$c] FAIL: empty sfs"; return 1; }
  "$REALSFS" saf2theta "$out.saf.idx" -sfs "$TH/baker.${c}.sfs" -outname "$TH/baker.${c}" -fold 1 2>>"$D/logs/saf.${c}.log"
  "$THETASTAT" do_stat "$TH/baker.${c}.thetas.idx" -outnames "$TH/baker.${c}.avg" 2>>"$D/logs/saf.${c}.log"
  "$THETASTAT" do_stat "$TH/baker.${c}.thetas.idx" -win 50000 -step 10000 -outnames "$TH/baker.${c}.win50k" 2>>"$D/logs/saf.${c}.log"
  [[ -s "$TH/baker.${c}.avg.pestPG" ]] && echo "[$(date)] [$c] DONE" || { echo "[$(date)] [$c] FAIL: no avg.pestPG"; return 1; }
}
export -f do_chrom
export ANGSD REALSFS THETASTAT REF BAM SITES SAF TH D THREADS
# all chroms except 2 (done); idempotent guard skips any already-complete
CHRS="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 Z"
echo "[$(date)] finalize sweep start (-P4 chroms, genome-wide .bin)"
printf '%s\n' $CHRS | xargs -P 4 -I{} bash -c 'do_chrom "$@"' _ {}
echo "[$(date)] finalize theta DONE"
