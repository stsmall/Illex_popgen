#!/usr/bin/env bash
# Incremental aggregation loop: re-aggregate all completed win50k.pestPG every ~12 min,
# keeping perpop_windows_expanded.tsv always valid. Exits after final aggregate once
# the driver marker appears.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/per_population
cd "$D"
PY=/home/ssmall/miniforge3/bin/python
LOG=logs/aggloop_expanded.log
echo "[$(date)] AGGLOOP START" > "$LOG"
while true; do
  echo "[$(date)] aggregate ($(ls thetas/*.win50k.pestPG 2>/dev/null | wc -l) win50k files)" >> "$LOG"
  "$PY" aggregate_theta_expanded.py >> "$LOG" 2>&1 || echo "[$(date)] aggregate error" >> "$LOG"
  [[ -f logs/EXPANDED_DRIVER_DONE ]] && break
  sleep 720
done
# final pass after driver done
echo "[$(date)] FINAL aggregate" >> "$LOG"
"$PY" aggregate_theta_expanded.py >> "$LOG" 2>&1
echo "[$(date)] AGGLOOP DONE" >> "$LOG"
touch logs/EXPANDED_AGG_DONE
