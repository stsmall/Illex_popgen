#!/usr/bin/env bash
# Per-chrom: SweepFinder2 unfolded CLR (thin -> -l, genome-wide background SFS) + RAiSD mu.
set -uo pipefail
c=$1
source "$(dirname "$0")/env.sh"
LEN=$(awk -v c="$c" '$1==c{print $2}' "$REF.fai")
G=$(( LEN / GRID_BP ));  [ "$G" -lt 1 ] && G=1        # SF2: 100kb grid
GR=$(( LEN / 10000 ));   [ "$GR" -lt 1 ] && GR=1      # RAiSD: 10kb grid
n0freq=$R/freq/$c.n$N0.freq
[ -s "$n0freq" ] || { echo "SCAN_FAIL $c (no projected freq)"; exit 1; }

# SweepFinder2: thin to 1 SNP/300bp (scan cost is O(SNP)/point), CLR vs genome-wide spectrum
"$PY" "$S/thin_freq.py" "$n0freq" "$R/freq/$c.thin.freq" bp 300 2>/dev/null
if "$SF2" -l "$G" "$R/freq/$c.thin.freq" "$R/sf2/genome.spect" "$R/sf2/$c.clr" > "$L/sf2_$c.log" 2>&1; then
  echo "SF2_OK $c ($(($(wc -l < "$R/sf2/$c.clr")-1)) pts)"
else echo "SF2_FAIL $c"; fi

# RAiSD composite mu on the polarized (REF=ancestral) VCF; -M 2 masks missing, -R gives factors
if ( cd "$R/raisd" && "$RAISD" -n "$c" -I "$R/polarized/$c.polarized.vcf.gz" -f -R -y 2 -M 2 -w 50 -G "$GR" \
      > "$L/raisd_$c.log" 2>&1 ); then
  echo "RAISD_OK $c ($(wc -l < "$R/raisd/RAiSD_Report.$c" 2>/dev/null || echo 0) pts)"
else echo "RAISD_FAIL $c"; fi
