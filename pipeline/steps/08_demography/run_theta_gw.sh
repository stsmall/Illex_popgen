#!/usr/bin/env bash
# Genome-wide diversity for P3 demography: baker-633 SAF -> realSFS -> theta (pi/thetaW/
# Tajima's D), per chromosome, FOLDED, accessible sites. Feeds windowed Tajima's D (3.1)
# + genome-wide SFS for Stairway (3.3). SHARED MACHINE: 4 chroms x -P6 ~= 24 cores (I/O-polite).
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
ANGSD=/home/ssmall/programs/angsd/angsd
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
THETASTAT=/home/ssmall/programs/angsd/misc/thetaStat
BBIN=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
D=$ANA/steps/08_demography
BAM=$ANA/steps/04_angsd_chr2/bamlists/bakerMeta_bamlist.txt   # 633 baker dedup
SAF=$D/saf; TH=$D/thetas; SITESD=$D/sites; mkdir -p "$SAF" "$TH" "$SITESD" "$D/logs"
touch "$REF.fai"   # ensure .fai newer than .fa (ANGSD is fatal otherwise)
THREADS=6

do_chrom() {
  local c=$1
  local sites="$SITESD/${c}.accessible.sites" out="$SAF/baker.${c}"
  [[ -s "$TH/baker.${c}.avg.pestPG" ]] && { echo "[$c] theta done"; return 0; }
  # per-chrom accessible sites (from accessible_sites.bed) + index
  if [[ ! -s "${sites}.idx" ]]; then
    awk -v C="$c" '$1==C{for(p=$2+1;p<=$3;p++) print C"\t"p}' "$ACCESS" > "$sites"
    "$ANGSD" sites index "$sites" >/dev/null 2>&1
  fi
  echo "[$(date)] [$c] SAF (633 baker, accessible, folded)"
  "$ANGSD" -b "$BAM" -sites "$sites" -r "$c" -anc "$REF" -ref "$REF" -out "$out" \
    -doSaf 1 -doCounts 1 -GL 1 -P "$THREADS" -minInd 1 -setMinDepthInd 1 \
    -minQ 20 -minMapQ 20 -remove_bads 1 -only_proper_pairs 1 > "$D/logs/saf.${c}.log" 2>&1
  "$REALSFS" "$out.saf.idx" -fold 1 -P "$THREADS" -tole 1e-8 -maxIter 200 > "$TH/baker.${c}.sfs" 2>>"$D/logs/saf.${c}.log"
  "$REALSFS" saf2theta "$out.saf.idx" -sfs "$TH/baker.${c}.sfs" -outname "$TH/baker.${c}" -fold 1 2>>"$D/logs/saf.${c}.log"
  "$THETASTAT" do_stat "$TH/baker.${c}.thetas.idx" -outnames "$TH/baker.${c}.avg" 2>>"$D/logs/saf.${c}.log"
  "$THETASTAT" do_stat "$TH/baker.${c}.thetas.idx" -win 50000 -step 10000 -outnames "$TH/baker.${c}.win50k" 2>>"$D/logs/saf.${c}.log"
  echo "[$(date)] [$c] DONE"
}
export -f do_chrom
export ANGSD REALSFS THETASTAT ACCESS REF BAM SAF TH SITESD D THREADS
# all chroms except 2 (bakerMeta chr2 SAF already exists in 04_angsd_chr2)
CHRS="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 Z"
echo "[$(date)] genome-wide baker SAF/theta start (-P4 chroms)"
printf '%s\n' $CHRS | xargs -P 4 -I{} bash -c 'do_chrom "$@"' _ {}
echo "[$(date)] genome-wide theta DONE"
