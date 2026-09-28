#!/usr/bin/env bash
# chr2 ReLERNN: AA homokaryotypes (background rho) + all-samples (inversion-suppressed), sequential on GPU2.
# Reuses the battle-tested per-tag driver run_one_tag.sh (UUID GPU pin, retry loops, --maskThresh 0.9).
set -uo pipefail
R=/sietch_colab/data_share/illex/popgen_data/analysis/steps/11_relernn
MASK=/sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/chr2_mask/chr2.mask.3state.bed
GPU=2
export SIM_EXTRA="--forceDiploid"   # chr2 subsets have samples fully-missing in some windows -> ReLERNN reads mixed ploidy (nSamps=-9); force diploid + md_mask handles missing
P=$R/chr2/relernn.progress
# chr2 inaccessible bed (ReLERNN -m expects INACCESSIBLE bases) + genome bed
awk 'BEGIN{OFS="\t"} $4=="inaccessible"{print $1,$2,$3}' "$MASK" > "$R/chr2/chr2.inacc.bed"
printf '2\t0\t119466599\n' > "$R/chr2/chr2.genome.bed"
setup(){ local tag=$1 vcf=$2
  mkdir -p "$R/$tag"
  ln -sf "$vcf" "$R/$tag/$tag.kept.vcf"
  cp -f "$R/chr2/chr2.genome.bed" "$R/$tag/genome.bed"
  cp -f "$R/chr2/chr2.inacc.bed"  "$R/$tag/inacc.bed"
}
setup chr2_AA  "$R/chr2/AA.chr2.vcf"
setup chr2_all "$R/chr2/all.chr2.vcf"
echo "CHR2_RELERNN_START $(date -u +%FT%TZ) gpu=$GPU" > "$P"
echo "AA_START $(date -u +%FT%TZ)" >> "$P"
bash "$R/run_one_tag.sh" chr2_AA "$GPU" >> "$P" 2>&1 && echo "AA_DONE $(date -u +%FT%TZ)" >> "$P" || echo "AA_FAIL $(date -u +%FT%TZ)" >> "$P"
echo "ALL_START $(date -u +%FT%TZ)" >> "$P"
bash "$R/run_one_tag.sh" chr2_all "$GPU" >> "$P" 2>&1 && echo "ALL_DONE $(date -u +%FT%TZ)" >> "$P" || echo "ALL_FAIL $(date -u +%FT%TZ)" >> "$P"
echo "CHR2_RELERNN_DONE $(date -u +%FT%TZ)" >> "$P"
