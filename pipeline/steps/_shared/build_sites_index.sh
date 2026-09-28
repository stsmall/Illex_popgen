#!/usr/bin/env bash
set -euo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
ANGSD=/home/ssmall/programs/angsd/angsd
SH=$ANA/steps/_shared
OUT=$SH/accessible_sites.sites
awk 'BEGIN{OFS="\t"} {for(p=$2+1;p<=$3;p++) print $1,p}' "$ACCESS" > "$OUT"
echo "[$(date)] expanded $(wc -l < "$OUT") positions; indexing..."
"$ANGSD" sites index "$OUT"
echo "[$(date)] DONE"; ls -la "$OUT"*
