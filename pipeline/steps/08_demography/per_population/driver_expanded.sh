#!/usr/bin/env bash
# Detached orchestrator: finish 4R gap, then compute full autosome set (excl chr2, chrZ).
# Strict core cap: CONC x P. Writes markers into logs/.
set -uo pipefail
D=/sietch_colab/data_share/illex/popgen_data/analysis/steps/08_demography/per_population
cd "$D"
LOG=logs/driver_expanded.log
echo "[$(date)] DRIVER START (load: $(cut -d' ' -f1 /proc/loadavg))" > "$LOG"
rm -f logs/EXPANDED_DRIVER_DONE

# Phase 0: finish the 4R gap (chroms 3,5,30). Reuse-aware; small batch.
echo "[$(date)] Phase0: 4R gap 3 5 30" >> "$LOG"
DIVS="4R" CHROMS="3 5 30" CONC=3 P=4 bash run_perpop_theta_expanded2.sh >> "$LOG" 2>&1

# Phase 1: full autosome set (interleaved sizes), CONC=4 x P=4 = 16 cores peak.
echo "[$(date)] Phase1: full autosome set" >> "$LOG"
CONC=4 P=4 bash run_perpop_theta_expanded2.sh >> "$LOG" 2>&1

echo "[$(date)] DRIVER DONE" >> "$LOG"
touch logs/EXPANDED_DRIVER_DONE
