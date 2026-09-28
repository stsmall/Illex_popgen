#!/usr/bin/env bash
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh

mkdir -p $WD/fvec/msc $WD/fvec/simstats $WD/fvec/train $WD/fvec/test

# --chrArmsForMasking requires comma-separated list (no spaces)
AUT_CSV="${AUT// /,}"

# Step 1: Split sims into hard/soft per subwindow position and write neutral msc
# (re-derives hard/soft from original RNG seeds; avoids any stdin-in-xargs issues)
echo "Splitting sweep sims by hard/soft type..."
$PY $WD/split_sweep_hardsoft.py

# Step 2: fvecSim for neutral
echo "Running fvecSim on neutral..."
$DSHIC fvecSim diploid $WD/fvec/msc/neutral.msc $WD/fvec/simstats/neutral.fvec \
  --totalPhysLen 44000 --numSubWins 11 \
  --maskFileName $WD/mask.fa --chrArmsForMasking $AUT_CSV \
  --unmaskedFracCutoff 0.25

# Step 3: fvecSim for hard and soft sweeps per subwindow
echo "Running fvecSim on hard/soft sweeps (11 subwindow positions each)..."
for s in $(seq 0 10); do
  echo "  fvecSim hard subwin $s..."
  $DSHIC fvecSim diploid $WD/fvec/msc/hard_${s}.msc $WD/fvec/simstats/hard_${s}.fvec \
    --totalPhysLen 44000 --numSubWins 11 \
    --maskFileName $WD/mask.fa --chrArmsForMasking $AUT_CSV \
    --unmaskedFracCutoff 0.25
  echo "  fvecSim soft subwin $s..."
  $DSHIC fvecSim diploid $WD/fvec/msc/soft_${s}.msc $WD/fvec/simstats/soft_${s}.fvec \
    --totalPhysLen 44000 --numSubWins 11 \
    --maskFileName $WD/mask.fa --chrArmsForMasking $AUT_CSV \
    --unmaskedFracCutoff 0.25
done

# Step 4: makeTrainingSets — separate hard/soft prefixes; sweep at subwin 5 (center)
# Files must end with _${i}.fvec; prefix is everything before _${i}.fvec
echo "Running makeTrainingSets..."
$DSHIC makeTrainingSets \
  $WD/fvec/simstats/neutral.fvec \
  $WD/fvec/simstats/soft \
  $WD/fvec/simstats/hard \
  5 \
  0,1,2,3,4,6,7,8,9,10 \
  $WD/fvec/train

echo "makeTrainingSets done. Files in $WD/fvec/train:"
ls -lh $WD/fvec/train/

# Step 5: 80/20 train/test split (20% → test)
echo "Creating 80/20 train/test split..."
$PY /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/split_traintest.py

echo "Train/test split done."
ls -lh $WD/fvec/train/ $WD/fvec/test/
