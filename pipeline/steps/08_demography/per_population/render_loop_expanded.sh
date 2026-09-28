#!/usr/bin/env bash
# Auto re-render the two diversity figures as the expanded compute extends.
# Re-renders whenever the number of fully-complete chromosomes (>=10 divisions)
# grows by >=3 since the last render, and once more after the driver finishes.
# Reads the table the aggloop keeps current; retries on transient read races.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/per_population
MAN=/sietch_colab/data_share/illex/popgen_data/analysis/manuscript
PY=/home/ssmall/miniforge3/envs/mkado-vcf/bin/python
export MPLCONFIGDIR=/dev/shm/mplcache
mkdir -p "$MPLCONFIGDIR"
cd "$MAN"
LOG=$D/logs/render_loop_expanded.log
echo "[$(date)] RENDER LOOP START" > "$LOG"

full_chroms() {
  for f in "$D"/thetas/*.win50k.pestPG; do b=$(basename "$f" .win50k.pestPG); echo "${b#*.}"; done \
    | sort | uniq -c | awk '$1>=10{n++} END{print n+0}'
}
render() {
  "$PY" plot_diversity.py >> "$LOG" 2>&1 && echo "[$(date)] rendered (full chroms=$1)" >> "$LOG" \
    || echo "[$(date)] render FAILED (retry next cycle)" >> "$LOG"
}

last=0
while true; do
  cur=$(full_chroms)
  if [[ $cur -ge $((last+3)) ]]; then render "$cur"; last=$cur; fi
  if [[ -f "$D/logs/EXPANDED_DRIVER_DONE" ]]; then
    echo "[$(date)] driver done -> final render (full chroms=$cur)" >> "$LOG"
    render "$cur"; break
  fi
  sleep 900
done
echo "[$(date)] RENDER LOOP DONE" >> "$LOG"
touch "$D/logs/EXPANDED_RENDER_DONE"
