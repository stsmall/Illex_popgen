#!/usr/bin/env bash
export CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=GPU-42ee9326-c26e-0977-b746-af7b9e4d57b5 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd /sietch_colab/data_share/illex/popgen_data/analysis/steps/09_momentsld/zA_males
/home/ssmall/programs/pg_gpu/.pixi/envs/default/bin/python compute_ld_zA.py Z > log_Z.txt 2>&1
/home/ssmall/programs/pg_gpu/.pixi/envs/default/bin/python compute_ld_zA.py A > log_A.txt 2>&1
touch LD_DONE
