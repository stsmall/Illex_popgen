#!/usr/bin/env bash
# Illex ensemble sweep scan on the BAKER SOURCE (full SFS): unfolded SweepFinder2 (CLR) +
# RAiSD (composite mu) + Bvalcalc (B-map).
# Run DETACHED:  setsid bash run_ensemble.sh > <progress> 2>&1 </dev/null &
set -uo pipefail
source "$(dirname "$0")/env.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PROG=$L/progress
echo "ENSEMBLE_START $(date -u +%FT%TZ)"

# ---- Bvalcalc genome-wide B-map (independent; background) ----
( if "$BVAL" --genome --params "$BVALDIR/Illex_Params.py" \
      --bedgff "$B/analysis/steps/16_bgspy/inputs/annot_cds.bed" \
      --chr_sizes "$BVALDIR/chr_sizes.csv" --rec_map "$BVALDIR/recmap.csv" \
      --gamma_dfe --out "$R/illex_Bmap.csv" --out_binsize 1000 --quiet > "$L/bvalcalc.log" 2>&1; then
    echo "BVALCALC_DONE $(date -u +%FT%TZ)" >> "$PROG"
  else echo "BVALCALC_FAIL" >> "$PROG"; fi ) &

# ---- Pass 1: polarize all chroms (baker ~18GB VCFs -> IO-bound, throttle) ----
echo "POLARIZE_START $(date -u +%FT%TZ)"
printf '%s\n' $CHROMS | xargs -P 8 -I{} bash "$S/polarize_chrom.sh" {}
echo "POLARIZE_DONE $(date -u +%FT%TZ)"

# ---- genome-wide background SFS (SF2 -f on pooled FULL projected n0 SNPs, unthinned) ----
first=$(echo $CHROMS | awk '{print $1}')
{ head -1 "$R/freq/$first.n$N0.freq"; for c in $CHROMS; do tail -n +2 "$R/freq/$c.n$N0.freq" 2>/dev/null; done; } > "$R/freq/genome.n$N0.freq"
"$SF2" -f "$R/freq/genome.n$N0.freq" "$R/sf2/genome.spect" > "$L/spect.log" 2>&1 \
  && echo "SPECT_DONE $(date -u +%FT%TZ)" || echo "SPECT_FAIL"

# ---- Pass 2: SF2 + RAiSD per chrom (parallel) ----
echo "SCAN_START $(date -u +%FT%TZ)"
printf '%s\n' $CHROMS | xargs -P 20 -I{} bash "$S/scan_chrom.sh" {}
echo "SCAN_DONE $(date -u +%FT%TZ)"

wait   # Bvalcalc background
echo "ENSEMBLE_DONE $(date -u +%FT%TZ)"
