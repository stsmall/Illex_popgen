#!/usr/bin/env bash
# Refit the windowed chr2 PCA (Supplementary Fig. S3) and the per-chromosome PCAs (eta^2 screen)
# on the 619 analysed individuals only (data/qc/karyotypes_619.tsv). The original fits used 653
# (windowed, chr2 callset) and 633 (per-chromosome, baker_633.txt) individuals, including the 12
# empty BAMs, the two contaminated individuals and, for chr2, 20 small_2026 samples.
#   windowed: bcftools subset of the chr2 callset, then pggpu_illex/popstats/scripts/winpca.py
#             unchanged (biallelic, MAF>0.05, LD prune, 200-kb windows, 100-kb step, k=10).
#   per-chromosome: thin_one.sh recipe (biallelic SNPs, MAF>0.05, +prune 1 SNP/10 kb) with the 619
#             list, then karyo_pca_619.py. 8 chromosomes at a time (shared machine).
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
M=$(cd "$HERE/../../.." && pwd)
ANA=/sietch_colab/data_share/illex/popgen_data/analysis
P=/sietch_colab/data_share/illex/popgen_data/pggpu_illex/popstats
BBIN=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin
PGPY=/home/ssmall/miniforge3/envs/varbuddy-pggpu/bin/python
O=${O:-$ANA/steps/05_karyo_pca_arms/refit_619}
GPU=${GPU:-1}
mkdir -p "$O/thinned" "$O/logs"
awk '!/^#/ && $1!="sample"{print $1}' "$M/data/qc/karyotypes_619.tsv" > "$O/keep_619.txt"
[ "$(wc -l < "$O/keep_619.txt")" -eq 619 ] || { echo "keep list is not 619"; exit 1; }
awk 'NR==FNR{k[$1];next} ($1 in k)' "$O/keep_619.txt" "$P/pops.tsv" > "$O/pops_619.tsv"

# ---- windowed chr2 PCA ----
if [ "${SKIP_WIN:-0}" != 1 ] && [ ! -s "$O/chr2_winpca_619_coords.npy" ]; then
  "$BBIN/bcftools" view -S "$O/keep_619.txt" -Oz -o "$O/chr2_619.vcf.gz" --threads 4 \
      "$ANA/steps/00_callset/filtered/2/variants_filt.vcf.gz" && "$BBIN/bcftools" index -t "$O/chr2_619.vcf.gz"
  "$BBIN/bcftools" query -l "$O/chr2_619.vcf.gz" > "$O/chr2_619_samples.txt"
  CUDA_VISIBLE_DEVICES=$GPU OMP_NUM_THREADS=4 "$PGPY" "$P/scripts/winpca.py" 2 "$O/chr2_619.vcf.gz" \
      "$O/pops_619.tsv" "$P/metadata.tsv" "$O/chr2_winpca_619" 200000 100000 > "$O/logs/winpca.log" 2>&1
  rm -f "$O/chr2_619.vcf.gz" "$O/chr2_619.vcf.gz.tbi"
fi

# ---- per-chromosome PCAs ----
thin_one() {
  local c=$1
  local out="$O/thinned/chr${c}.thin.vcf.gz"
  [[ -s "$out.tbi" ]] && return 0
  "$BBIN/bcftools" view -S "$O/keep_619.txt" -m2 -M2 -v snps -q 0.05:minor \
      "$ANA/steps/00_callset/filtered/${c}/variants_filt.vcf.gz" -Ou 2>/dev/null \
    | "$BBIN/bcftools" +prune -w 10kb -n 1 -Oz -o "$out" 2>/dev/null
  "$BBIN/bcftools" index -t "$out" && echo "chr$c: $("$BBIN/bcftools" index -n "$out") SNPs"
}
export -f thin_one; export BBIN O ANA
{ seq 1 45; echo Z; } | xargs -P 8 -I{} bash -c 'thin_one "$@"' _ {} > "$O/logs/thin.log" 2>&1
CUDA_VISIBLE_DEVICES=$GPU OMP_NUM_THREADS=4 "$PGPY" "$HERE/karyo_pca_619.py" "$O" > "$O/logs/karyo_pca.log" 2>&1
[ -s "$O/karyo_pca_coords_619.tsv" ] && echo DONE > "$O/ALL_DONE"
