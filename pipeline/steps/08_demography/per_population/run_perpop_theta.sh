#!/usr/bin/env bash
# Per-(division x chromosome) ANGSD SAF -> realSFS -> thetas (pi/thetaW/Tajima's D),
# FOLDED, accessible sites. Mirrors 08_demography/run_theta_gw2.sh exactly, split by
# (division,chrom) for parallelism. Windowed do_stat (50kb non-overlap) -> per-window
# pi/TajD for box plots + per-chrom facet. SHARED MACHINE: CONC*P threads for my jobs.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
ANGSD=/home/ssmall/programs/angsd/angsd
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
THETASTAT=/home/ssmall/programs/angsd/misc/thetaStat
D=$ANA/steps/08_demography/per_population
SITES=$ANA/steps/_shared/accessible_sites.sites   # genome-wide, indexed (.bin/.idx)
SAF=$D/saf; TH=$D/thetas; mkdir -p "$SAF" "$TH" "$D/logs"
touch "$REF.fai"
P=${P:-4}

do_job() {
  local job=$1; local d=${job%%:*}; local c=${job##*:}
  local bam="$D/bamlists_sub/${d}.bamlist.txt" out="$SAF/${d}.${c}"
  [[ -s "$TH/${d}.${c}.avg.pestPG" && -s "$TH/${d}.${c}.win50k.pestPG" ]] && { echo "[$d $c] done, skip"; return 0; }
  "$ANGSD" -b "$bam" -sites "$SITES" -r "$c" -anc "$REF" -ref "$REF" -out "$out" \
    -doSaf 1 -doCounts 1 -GL 1 -P "$P" -minInd 1 -setMinDepthInd 1 \
    -minQ 20 -minMapQ 20 -remove_bads 1 -only_proper_pairs 1 > "$D/logs/${d}.${c}.log" 2>&1
  [[ -s "$out.saf.idx" ]] || { echo "[$d $c] FAIL saf.idx"; return 1; }
  "$REALSFS" "$out.saf.idx" -fold 1 -P "$P" -tole 1e-8 -maxIter 200 > "$TH/${d}.${c}.sfs" 2>>"$D/logs/${d}.${c}.log"
  [[ -s "$TH/${d}.${c}.sfs" ]] || { echo "[$d $c] FAIL sfs"; return 1; }
  "$REALSFS" saf2theta "$out.saf.idx" -sfs "$TH/${d}.${c}.sfs" -outname "$TH/${d}.${c}" -fold 1 2>>"$D/logs/${d}.${c}.log"
  "$THETASTAT" do_stat "$TH/${d}.${c}.thetas.idx" -outnames "$TH/${d}.${c}.avg" 2>>"$D/logs/${d}.${c}.log"
  "$THETASTAT" do_stat "$TH/${d}.${c}.thetas.idx" -win 50000 -step 50000 -type 0 -outnames "$TH/${d}.${c}.win50k" 2>>"$D/logs/${d}.${c}.log"
  # prune bulky per-job SAF to save disk (thetas kept)
  rm -f "$out.saf.gz" "$out.saf.pos.gz" "$out.saf.idx"
  [[ -s "$TH/${d}.${c}.avg.pestPG" ]] && echo "[$(date +%H:%M)] [$d $c] DONE" || { echo "[$d $c] FAIL avg"; return 1; }
}
export -f do_job
export ANGSD REALSFS THETASTAT REF SITES SAF TH D P
DIVS="${DIVS:-6C 6B 6A 4X 4W 3O 4T 4R 4S 3K}"
CHROMS="${CHROMS:-1 3 5 9 13 21}"
CONC="${CONC:-12}"
JOBS=""
for d in $DIVS; do for c in $CHROMS; do JOBS="$JOBS $d:$c"; done; done
echo "[$(date)] perpop theta start: $(echo $JOBS|wc -w) jobs, CONC=$CONC, P=$P"
printf '%s\n' $JOBS | xargs -P "$CONC" -I{} bash -c 'do_job "$@"' _ {}
echo "[$(date)] perpop theta ALL DONE"
