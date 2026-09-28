#!/usr/bin/env bash
# Z maps on GPU1. GPU0 lost sole-tenancy (jaxeng JAX kernel co-resident -> cuDNN CudnnRNNV3
# NOT_SUPPORTED on GRU train). GPU1 is nearly free (only a 416MB jaxeng device handle).
# run_one_tag.sh passes --gpuID 1 + UUID pinning + CUDA_DEVICE_ORDER=PCI_BUS_ID (ordering-immune).
# Sims for both Z tags are already cached -> goes straight to TRAIN (reuse existing sims).
set -uo pipefail
RD=/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn
echo "[$(date)] Z-on-GPU1 driver start"
bash $RD/run_one_tag.sh Z_male 1
bash $RD/run_one_tag.sh Z_female 1
echo "[$(date)] Z-on-GPU1 driver DONE"
