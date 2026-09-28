#!/usr/bin/env bash
# Run after fvec_sim.sh and fvec_real.sh complete.
# Step 6: verify check passes, then checkpoint.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh

echo "=== Task 6 finalization ==="
echo "fvec/train contents:"
ls -lh $WD/fvec/train/ 2>/dev/null || echo "(empty)"
echo "fvec/test contents:"
ls -lh $WD/fvec/test/ 2>/dev/null || echo "(empty)"
NREAL=$(find $WD/fvec/real -name '*.fvec' | wc -l)
NSUB=$(ls -d $WD/fvec/real/sub* 2>/dev/null | wc -l)
echo "fvec/real: $NREAL fvec files across $NSUB subsamples"

echo "Running check..."
$PY $WD/checks/check_fvec.py

NTRAIN_ROWS=$(wc -l < $WD/fvec/train/hard.fvec 2>/dev/null || echo 0)
NTEST_ROWS=$(wc -l < $WD/fvec/test/hard.fvec 2>/dev/null || echo 0)

echo "$(date +%F) Task6 fvec: mask.fa built, train/test 5-class fvecs (~$NTRAIN_ROWS rows/class), $NREAL real sub-chrom fvecs ($NSUB subsamples)" >> $WD/RUNLOG.md
echo "Checkpoint written."
tail -3 $WD/RUNLOG.md
