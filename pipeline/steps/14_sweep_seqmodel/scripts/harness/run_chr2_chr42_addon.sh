#!/usr/bin/env bash
# Add-on ensemble scan (RAiSD + SF2) for chr2 + chr42, which the original run excluded (stale
# CHROMS list). Reuses the existing genome-wide background spectrum (results/sf2/genome.spect) —
# 2 more chroms negligibly change the genome SFS. polarize -> scan, per chrom, 2-way parallel.
# NOTE chr2 = inversion chromosome: its sweep calls are inversion-confounded; interpret collinear
# regions only (the inversion body's signal is balancing selection, §3, not classic sweeps).
# Detached:  setsid bash run_chr2_chr42_addon.sh >/dev/null 2>&1 </dev/null &
set -uo pipefail
ES=/sietch_colab/data_share/illex/popgen_data/analysis/steps/17_ensemble_scan
S=$ES/scripts
mkdir -p "$ES"/results/{polarized,freq,sf2,raisd} "$ES"/logs
PROG=$ES/logs/chr2_42_addon.progress; : > "$PROG"
say(){ echo "$(date -u +%FT%TZ) $*" >> "$PROG"; }
export S PROG
say "ADDON_START chr2 chr42 (reusing genome.spect)"
do_one(){
  local c=$1
  echo "$(date -u +%FT%TZ) POLARIZE $c" >> "$PROG"
  if bash "$S/polarize_chrom.sh" "$c" >> "$PROG" 2>&1; then
    echo "$(date -u +%FT%TZ) SCAN $c" >> "$PROG"
    bash "$S/scan_chrom.sh" "$c" >> "$PROG" 2>&1
  fi
  echo "$(date -u +%FT%TZ) CHROM_DONE $c (raisd=$([ -s "$ES/results/raisd/RAiSD_Report.$c" ]&&echo Y||echo N) sf2=$([ -s "$ES/results/sf2/$c.clr" ]&&echo Y||echo N))" >> "$PROG"
}
export -f do_one; export ES
printf '%s\n' 2 42 | xargs -P 2 -I{} bash -c 'do_one "$@"' _ {}
say "ADDON_DONE"
