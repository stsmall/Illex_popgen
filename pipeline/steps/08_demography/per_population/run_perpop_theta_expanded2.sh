#!/usr/bin/env bash
# REUSE-AWARE expanded per-(division x chromosome) ANGSD SAF -> realSFS -> thetas.
# IDENTICAL method/params to run_perpop_theta.sh / run_perpop_theta_expanded.sh.
# Improvements: reuses existing SAF (.saf.idx) and SFS (.sfs, non-empty) intermediates
# instead of recomputing; only cleans SAF files after the job fully completes.
# Strict core cap: CONC concurrent jobs x P threads (default 4x4 = 16 cores).
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
ANGSD=/home/ssmall/programs/angsd/angsd
REALSFS=/home/ssmall/programs/angsd/misc/realSFS
THETASTAT=/home/ssmall/programs/angsd/misc/thetaStat
D=$ANA/steps/08_demography/per_population
SITES=$ANA/steps/_shared/accessible_sites.sites
SAF=$D/saf; TH=$D/thetas; mkdir -p "$SAF" "$TH" "$D/logs"
touch "$REF.fai"
P=${P:-4}

do_job() {
  local job=$1; local d=${job%%:*}; local c=${job##*:}
  local bam="$D/bamlists_sub/${d}.bamlist.txt" out="$SAF/${d}.${c}"
  local log="$D/logs/${d}.${c}.log"
  # already fully done?
  [[ -s "$TH/${d}.${c}.avg.pestPG" && -s "$TH/${d}.${c}.win50k.pestPG" ]] && { echo "[$d $c] done, skip"; return 0; }
  # reuse SAF if present, else compute
  if [[ ! -s "$out.saf.idx" ]]; then
    "$ANGSD" -b "$bam" -sites "$SITES" -r "$c" -anc "$REF" -ref "$REF" -out "$out" \
      -doSaf 1 -doCounts 1 -GL 1 -P "$P" -minInd 1 -setMinDepthInd 1 \
      -minQ 20 -minMapQ 20 -remove_bads 1 -only_proper_pairs 1 > "$log" 2>&1
    [[ -s "$out.saf.idx" ]] || { echo "[$d $c] FAIL saf.idx"; return 1; }
  else
    echo "[$d $c] reuse SAF" >> "$log"
  fi
  # reuse SFS if non-empty, else compute
  if [[ ! -s "$TH/${d}.${c}.sfs" ]]; then
    "$REALSFS" "$out.saf.idx" -fold 1 -P "$P" -tole 1e-8 -maxIter 200 > "$TH/${d}.${c}.sfs" 2>>"$log"
    [[ -s "$TH/${d}.${c}.sfs" ]] || { echo "[$d $c] FAIL sfs"; return 1; }
  else
    echo "[$d $c] reuse SFS" >> "$log"
  fi
  "$REALSFS" saf2theta "$out.saf.idx" -sfs "$TH/${d}.${c}.sfs" -outname "$TH/${d}.${c}" -fold 1 2>>"$log"
  "$THETASTAT" do_stat "$TH/${d}.${c}.thetas.idx" -outnames "$TH/${d}.${c}.avg" 2>>"$log"
  "$THETASTAT" do_stat "$TH/${d}.${c}.thetas.idx" -win 50000 -step 50000 -type 0 -outnames "$TH/${d}.${c}.win50k" 2>>"$log"
  if [[ -s "$TH/${d}.${c}.avg.pestPG" && -s "$TH/${d}.${c}.win50k.pestPG" ]]; then
    rm -f "$out.saf.gz" "$out.saf.pos.gz" "$out.saf.idx"
    echo "[$(date +%H:%M)] [$d $c] DONE"
  else
    echo "[$d $c] FAIL avg/win"; return 1
  fi
}
export -f do_job
export ANGSD REALSFS THETASTAT REF SITES SAF TH D P

DIVS="${DIVS:-3K 3O 4R 4S 4T 4W 4X 6A 6B 6C}"
# CHROMS order chosen to interleave size classes so incremental snapshots span sizes.
CHROMS="${CHROMS:-6 40 25 11 41 34 9 44 16 4 22 27 14 39 20 10 35 24 12 36 43 7 37 18 8 38 29 31 28 15 23 32 19 26 42}"
CONC="${CONC:-4}"
# Build job list chrom-major so whole chroms finish together (good for incremental aggregation).
JOBS=""
for c in $CHROMS; do for d in $DIVS; do JOBS="$JOBS $c:$d"; done; done
# NOTE: job token is c:d here; do_job expects d:c, so swap.
JOBS2=""
for j in $JOBS; do c=${j%%:*}; d=${j##*:}; JOBS2="$JOBS2 $d:$c"; done
echo "[$(date)] REUSE expanded perpop theta start: $(echo $JOBS2|wc -w) jobs, CONC=$CONC, P=$P (peak cores=$((CONC*P)))"
printf '%s\n' $JOBS2 | xargs -P "$CONC" -I{} bash -c 'do_job "$@"' _ {}
echo "[$(date)] REUSE expanded perpop theta ALL DONE"
