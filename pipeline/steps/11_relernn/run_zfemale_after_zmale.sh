#!/usr/bin/env bash
# Run Z_female on GPU0 AFTER Z_male completes (BSCORRECTED appears). GPU0 is the reliable
# sole-tenant-ish card; GPU1 flickers constantly (adkern) and tripped cuDNN 3x for Z_female.
# run_one_tag.sh Z_female 0 has 8x TRAIN auto-retry to ride out any transient GPU0 flicker.
set -uo pipefail
RD=/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn
echo "[$(date)] waiting for Z_male BSCORRECTED before starting Z_female on GPU0..."
until find $RD/run_Z_male/proj -name '*BSCORRECTED.txt' 2>/dev/null | grep -q .; do sleep 300; done
echo "[$(date)] Z_male done -> starting Z_female on GPU0"
bash $RD/run_one_tag.sh Z_female 0
echo "[$(date)] Z_female (GPU0) driver finished"
