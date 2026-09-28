#!/usr/bin/env bash
# pilot_status.sh -- one-shot status of the soft pilot campaign, for the cron monitor.
# Prints machine-readable facts; the caller judges DONE / STUCK / RUNNING / INCOMPLETE.
OUT=/sietch_colab/data_share/illex/popgen_data/analysis/steps/14_sweep_seqmodel/results/pilot/train_ms
LOG=$OUT/campaign_soft.log
# run_campaign soft alive? generate.py soft sims running? (comm-filtered -> no self-match)
alive=$(ps -u ssmall -o comm=,args= 2>/dev/null | awk '$1 ~ /python/ && /run_campaign\.py/ && /--klass soft/{n++} END{print (n>0)?"alive":"gone"}')
nsim=$(ps -u ssmall -o comm=,args= 2>/dev/null | awk '$1 ~ /python/ && /generate\.py/ && /--klass soft/{n++} END{print n+0}')
ncomb=$(ls "$OUT"/soft_*.msOut.gz 2>/dev/null | wc -l)
echo "CAMPAIGN_PROC: $alive"
echo "SOFT_SIMS_RUNNING: $nsim"
echo "SOFT_COMBINED: $ncomb/11"
echo "SOFT_SIM_DIRS: $(ls -d "$OUT"/sims/soft/s* 2>/dev/null | wc -l)/40 done-or-running"
if [ -f "$LOG" ]; then
  echo "LOG_MTIME_AGO_MIN: $(( ( $(date +%s) - $(stat -c %Y "$LOG") ) / 60 ))"
  echo "FAILURES_LINES: $(grep -c 'FAILURES' "$LOG" 2>/dev/null)"
  echo "LOG_TAIL:"; tail -4 "$LOG" 2>/dev/null | sed 's/^/  /'
else
  echo "LOG: MISSING (campaign not started or died before first write)"
fi