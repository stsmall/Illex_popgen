#!/usr/bin/env bash
# Per-division ANGSD SAF -> realSFS -> thetas (pi/thetaW/TajD), genome-wide, FOLDED,
# accessible sites. Mirrors 08_demography/run_theta_gw2.sh but per NAFO division.
# Windowed do_stat (win50k/step10k) gives per-window pi/TajD (with Chr col) for
# box plots (main) + per-chrom facet (supp). SHARED MACHINE: keep footprint modest.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
ANGSD=/home/ssmall/programs/angsd/angsd
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
THETASTAT=/home/ssmall/programs/angsd/misc/thetaStat
D=$ANA/steps/08_demography/per_population
SITES=$ANA/steps/_shared/accessible_sites.sites   # genome-wide, indexed (.bin/.idx)
SAF=$D/saf; TH=$D/thetas; mkdir -p "$SAF" "$TH" "$D/logs"
touch "$REF.fai"
THREADS=${THREADS:-4}

do_div() {
  local d=$1
  local bam="$D/bamlists/${d}.bamlist.txt" out="$SAF/${d}"
  [[ -s "$TH/${d}.avg.pestPG" ]] && { echo "[$d] theta done, skip"; return 0; }
  echo "[$(date)] [$d] SAF start (N=$(wc -l < "$bam"))"
  "$ANGSD" -b "$bam" -sites "$SITES" -anc "$REF" -ref "$REF" -out "$out" \
    -doSaf 1 -doCounts 1 -GL 1 -P "$THREADS" -minInd 1 -setMinDepthInd 1 \
    -minQ 20 -minMapQ 20 -remove_bads 1 -only_proper_pairs 1 > "$D/logs/saf.${d}.log" 2>&1
  [[ -s "$out.saf.idx" ]] || { echo "[$(date)] [$d] FAIL: empty saf.idx"; return 1; }
  "$REALSFS" "$out.saf.idx" -fold 1 -P "$THREADS" -tole 1e-8 -maxIter 200 > "$TH/${d}.sfs" 2>>"$D/logs/saf.${d}.log"
  [[ -s "$TH/${d}.sfs" ]] || { echo "[$(date)] [$d] FAIL: empty sfs"; return 1; }
  "$REALSFS" saf2theta "$out.saf.idx" -sfs "$TH/${d}.sfs" -outname "$TH/${d}" -fold 1 2>>"$D/logs/saf.${d}.log"
  "$THETASTAT" do_stat "$TH/${d}.thetas.idx" -outnames "$TH/${d}.avg" 2>>"$D/logs/saf.${d}.log"
  "$THETASTAT" do_stat "$TH/${d}.thetas.idx" -win 50000 -step 50000 -type 0 -outnames "$TH/${d}.win50k" 2>>"$D/logs/saf.${d}.log"
  [[ -s "$TH/${d}.avg.pestPG" ]] && echo "[$(date)] [$d] DONE" || { echo "[$(date)] [$d] FAIL: no avg.pestPG"; return 1; }
}
export -f do_div
export ANGSD REALSFS THETASTAT REF SITES SAF TH D THREADS
DIVS="${DIVS:-6C 6B 6A 4X 4W 3O 4T 4R 4S 3K}"
CONC="${CONC:-3}"
echo "[$(date)] per-division SAF/theta start (conc=$CONC, P=$THREADS): $DIVS"
printf '%s\n' $DIVS | xargs -P "$CONC" -I{} bash -c 'do_div "$@"' _ {}
echo "[$(date)] per-division theta DONE"
