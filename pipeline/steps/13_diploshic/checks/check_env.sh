#!/usr/bin/env bash
set -euo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh
$DISCOAL 5 1 100 -t 5 >/dev/null
echo "discoal OK"
$DSHIC --help >/dev/null 2>&1
echo "diploSHIC OK"
$PY -c "import allel, tensorflow, numpy, pandas"
echo "deps OK"
echo "ALL ENV CHECKS PASSED"
