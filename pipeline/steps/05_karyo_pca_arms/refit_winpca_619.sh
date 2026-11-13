#!/usr/bin/env bash
# Windowed chr2 PCA (Supplementary Fig. S1) on the 619 analysed individuals. Same as the windowed step of
# refit_pca_619.sh, but run in the pg_gpu pixi environment: pg_gpu's local_pca needs cusolver (eigh), which the
# varbuddy-pggpu conda env lacks (ImportError libcusolver.so.11). The original genome-wide scan
# (pggpu_illex/popstats/scripts/run_winpca_all.sh) also ran winpca.py through pixi.
set -uo pipefail
ANA=/sietch_colab/data_share/illex/popgen_data/analysis
P=/sietch_colab/data_share/illex/popgen_data/pggpu_illex/popstats
BBIN=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
PGREPO=/home/ssmall/programs/pg_gpu
O=${O:-$ANA/steps/05_karyo_pca_arms/refit_619}
GPU=${GPU:-1}
[ "$(wc -l < "$O/keep_619.txt")" -eq 619 ] || { echo "run refit_pca_619.sh first (keep list)"; exit 1; }
if [ ! -s "$O/chr2_619.vcf.gz.tbi" ]; then
  "$BBIN/bcftools" view -S "$O/keep_619.txt" -Oz -o "$O/chr2_619.vcf.gz" --threads 4 \
      "$ANA/steps/00_callset/filtered/2/variants_filt.vcf.gz" && "$BBIN/bcftools" index -t "$O/chr2_619.vcf.gz"
fi
"$BBIN/bcftools" query -l "$O/chr2_619.vcf.gz" > "$O/chr2_619_samples.txt"
( cd "$PGREPO" && CUDA_VISIBLE_DEVICES=$GPU OMP_NUM_THREADS=4 pixi run python "$P/scripts/winpca.py" 2 \
    "$O/chr2_619.vcf.gz" "$O/pops_619.tsv" "$P/metadata.tsv" "$O/chr2_winpca_619" 200000 100000 ) \
    > "$O/logs/winpca.log" 2>&1
if [ -s "$O/chr2_winpca_619_coords.npy" ]; then
  rm -f "$O/chr2_619.vcf.gz" "$O/chr2_619.vcf.gz.tbi"; echo DONE > "$O/WIN_DONE"
else
  echo FAILED > "$O/WIN_FAILED"
fi
