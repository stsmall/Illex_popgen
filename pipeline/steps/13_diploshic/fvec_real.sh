#!/usr/bin/env bash
set -euo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh

K=50
NDIP=350

# Step 1: Generate fixed subsample lists (seeded, reproducible, without replacement)
echo "Generating $K subsample lists of $NDIP samples..."
$PY $WD/gen_subsamples.py $K $NDIP

# Step 2: Verify samples file exists and has enough samples
NTOTAL=$(wc -l < $WD/all_samples.txt)
echo "Total samples: $NTOTAL; drawing $NDIP without replacement per subsample"
if [ $NTOTAL -lt $NDIP ]; then
  echo "ERROR: not enough samples ($NTOTAL < $NDIP)" >&2
  exit 1
fi

mkdir -p $WD/fvec/real

# Step 3: K*43 fvecVcf jobs at P16 (shared box — do not raise)
# Use external fvec_one.sh to avoid stdin-in-xargs bash-function gotcha
echo "Running fvecVcf: $K subsamples x 43 chromosomes = $((K * 43)) jobs at P16..."
chmod +x $WD/fvec_one.sh

for k in $(seq 0 $((K-1))); do
  for c in $AUT; do
    echo "$k $c"
  done
done | xargs -P16 -n2 bash $WD/fvec_one.sh

NFVEC=$(find $WD/fvec/real -name '*.fvec' | wc -l)
echo "Done: $NFVEC subsample-chrom fvecs produced (expect ~$((K * 43)) = $((K * 43)))"
