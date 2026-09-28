#!/usr/bin/env bash
# GONE2 chromosome-block bootstrap. For each replicate: build a resampled VCF
# (chromosomes drawn with replacement, relabelled to fresh linkage groups), run
# GONE2 with the same core flags as the point estimate (-g 3 -r 1) plus -s 50000
# (random 50k-SNP subsample; verified to reproduce the full-SNP point-estimate
# trajectory essentially exactly -- ratio 1.00 at every generation, t201 test --
# while cutting runtime from ~45 min to ~3.5 min/replicate). Save each *_GONE2_Ne
# trajectory, then delete the big replicate VCF. 4 parallel x 24 threads = 96 cores
# (under the 100-core shared-machine cap).
set -uo pipefail
BD=/sietch_colab/data_share/illex/popgen_data/analysis/steps/10_gone2/bootstrap
GONE2=/home/ssmall/programs/GONE2/gone2
PY=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/python
NREP=${1:-30}
THREADS=24
PAR=4
SNPS=50000
mkdir -p "$BD/reps" "$BD/ne" "$BD/logs"

one() {
  local rep=$1
  local vcf="$BD/reps/rep${rep}.vcf"
  local log="$BD/logs/rep${rep}.log"
  local ne="$BD/ne/rep${rep}_GONE2_Ne"
  [[ -s "$ne" ]] && { echo "[rep $rep] already done, skip"; return 0; }
  echo "[$(date +%H:%M:%S)] [rep $rep] build vcf"
  "$PY" "$BD/make_boot_vcf.py" "$rep" "$vcf" >>"$log" 2>&1 || { echo "[rep $rep] BUILD FAIL"; return 1; }
  echo "[$(date +%H:%M:%S)] [rep $rep] run gone2"
  "$GONE2" -g 3 -r 1 -s "$SNPS" -t "$THREADS" "$vcf" >>"$log" 2>&1
  if [[ -s "${vcf}_GONE2_Ne" ]]; then
    cp "${vcf}_GONE2_Ne" "$ne"
    echo "[$(date +%H:%M:%S)] [rep $rep] DONE -> $ne"
  else
    echo "[$(date +%H:%M:%S)] [rep $rep] GONE2 FAIL (no _GONE2_Ne)"
  fi
  # cleanup big files, keep the log
  rm -f "$vcf" "${vcf}"_GONE2_* "${vcf}"_GONE_* 2>/dev/null
}
export -f one
export BD GONE2 PY THREADS SNPS

echo "[$(date)] GONE2 bootstrap: $NREP reps, $PAR parallel x $THREADS threads, -s $SNPS"
seq 1 "$NREP" | xargs -P "$PAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date)] ALL DONE. Ne trajectories in $BD/ne/"
ls "$BD/ne/" | wc -l
