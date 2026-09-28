#!/usr/bin/env bash
# One-shot status of the ensemble scan, for the cron monitor.
source "$(dirname "$0")/env.sh"
alive=$(ps -u ssmall -o comm=,args= 2>/dev/null | awk '$1=="bash" && /run_ensemble/{n++} END{print (n>0)?"alive":"gone"}')
echo "RUN_ENSEMBLE: $alive"
echo "FREQ: $(ls "$R"/freq/*.sf2.freq 2>/dev/null | grep -c . )/43"
echo "SF2_CLR: $(ls "$R"/sf2/*.clr 2>/dev/null | grep -c . )/43"
echo "RAISD: $(ls "$R"/raisd/RAiSD_Report.* 2>/dev/null | grep -c . )/43"
echo "BMAP: $([ -s "$R/illex_Bmap.csv" ] && echo done || echo pending)"
echo "FAILS: $(grep -cE '_FAIL' "$L/progress" 2>/dev/null)"
echo "SENTINELS: $(grep -oE 'POLARIZE_DONE|SPECT_DONE|SCAN_DONE|ENSEMBLE_DONE|BVALCALC_DONE' "$L/progress" 2>/dev/null | tr '\n' ' ')"
echo "TAIL:"; grep -E 'START|_DONE|_FAIL' "$L/progress" 2>/dev/null | tail -5 | sed 's/^/  /'