#!/usr/bin/env bash
# Genome-wide diversity for I. argentinus (small_2026, 10 low-cov dedup BAMs mapped to
# Illex_F24 illecebrosus reference). EXACT replicate of the illecebrosus route-A pipeline
# (08_demography/run_theta_gw2.sh): SAF -> realSFS (folded) -> saf2theta -> thetaStat,
# per chromosome, FOLDED, restricted to the genome-wide accessible-sites index.
# Same ANGSD flags: -doSaf 1 -doCounts 1 -GL 1 -minInd 1 -setMinDepthInd 1 -minQ 20
# -minMapQ 20 -remove_bads 1 -only_proper_pairs 1, anc=ref (folded).
# SHARED MACHINE: 8 chroms x -P6 ~= 48 threads.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
ANGSD=/home/ssmall/miniforge3/envs/angsd/bin/angsd
REALSFS=/home/ssmall/miniforge3/envs/angsd/bin/realSFS
THETASTAT=/home/ssmall/miniforge3/envs/angsd/bin/thetaStat
D=$ANA/steps/08_demography/argentinus
BAM=$D/argentinus_bamlist.txt                                  # 10 argentinus dedup
SITES=$ANA/steps/_shared/accessible_sites.sites               # genome-wide, indexed (.bin/.idx)
SAF=$D/saf; TH=$D/thetas; mkdir -p "$SAF" "$TH" "$D/logs"
touch "$REF.fai"
THREADS=6

do_chrom() {
  local c=$1
  local out="$SAF/arg.${c}"
  [[ -s "$TH/arg.${c}.avg.pestPG" ]] && { echo "[$c] theta done, skip"; return 0; }
  echo "[$(date)] [$c] SAF (10 argentinus, genome-wide .bin sites, folded)"
  "$ANGSD" -b "$BAM" -sites "$SITES" -r "$c" -anc "$REF" -ref "$REF" -out "$out" \
    -doSaf 1 -doCounts 1 -GL 1 -P "$THREADS" -minInd 1 -setMinDepthInd 1 \
    -minQ 20 -minMapQ 20 -remove_bads 1 -only_proper_pairs 1 > "$D/logs/saf.${c}.log" 2>&1
  [[ -s "$out.saf.idx" ]] || { echo "[$(date)] [$c] FAIL: empty saf.idx"; return 1; }
  "$REALSFS" "$out.saf.idx" -fold 1 -P "$THREADS" -tole 1e-8 -maxIter 200 > "$TH/arg.${c}.sfs" 2>>"$D/logs/saf.${c}.log"
  [[ -s "$TH/arg.${c}.sfs" ]] || { echo "[$(date)] [$c] FAIL: empty sfs"; return 1; }
  "$REALSFS" saf2theta "$out.saf.idx" -sfs "$TH/arg.${c}.sfs" -outname "$TH/arg.${c}" -fold 1 2>>"$D/logs/saf.${c}.log"
  "$THETASTAT" do_stat "$TH/arg.${c}.thetas.idx" -outnames "$TH/arg.${c}.avg" 2>>"$D/logs/saf.${c}.log"
  [[ -s "$TH/arg.${c}.avg.pestPG" ]] && echo "[$(date)] [$c] DONE" || { echo "[$(date)] [$c] FAIL: no avg.pestPG"; return 1; }
}
export -f do_chrom
export ANGSD REALSFS THETASTAT REF BAM SITES SAF TH D THREADS
# all 46 chroms; pooling (excl chr2 inversion, chr42 sex, chrZ) done in summarize step
CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 Z"
echo "[$(date)] argentinus genome-wide SAF/theta start (-P8 chroms x -P6)"
printf '%s\n' $CHRS | xargs -P 8 -I{} bash -c 'do_chrom "$@"' _ {}
echo "[$(date)] argentinus theta DONE"
