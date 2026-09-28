#!/usr/bin/env bash
# Parallel merge driver. SHARED MACHINE: hard cap 100 cores. Here: 12 concurrent
# chroms x 2 compression threads each ~= ~24-36 cores. Do NOT raise -P beyond
# what keeps (concurrency x threads) well under 100.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
cd "$ANA/steps/00_callset"
CONCURRENCY=12
echo "[$(date)] parallel merge start (-P $CONCURRENCY, 2 threads/job)"
printf '%s\n' $CHR_ALL | xargs -P "$CONCURRENCY" -I{} bash merge_one.sh {}
echo "[$(date)] parallel merge driver done"
