#!/usr/bin/env bash
# Prepare mkado inputs: subset each per-chrom VCF to CDS regions and strip FORMAT to GT only,
# then concat into genome-wide CDS-only bgzipped+tabix VCFs. mkado only uses CDS sites + GT,
# so this is correct and ~1000x smaller than the full 1.1TB callset.
set -euo pipefail

ROOT=/sietch_colab/data_share/illex/popgen_data/mkado_illex
BAKER=/sietch_colab/data_share/illex/popgen_data/seq_data/baker_2025/vcfs
ARG=/sietch_colab/data_share/illex/popgen_data/seq_data/small_2026/vcfs/arg
IN=$ROOT/inputs
CDS=$ROOT/cds
BED=$IN/cds.bed
mkdir -p "$IN" "$CDS"

CHROMS=$(for i in $(seq 1 45); do echo $i; done; echo Z)

# subset one chrom: view CDS regions (index random-access) + keep only GT, drop INFO
subset_chrom() {
  local src=$1 out=$2 bed=$3
  bcftools view -R "$bed" -Ou "$src" \
    | bcftools annotate -x 'INFO,^FORMAT/GT' -Oz -o "$out"
  tabix -f -p vcf "$out"
}
export -f subset_chrom

echo "[$(date)] === Ingroup baker: per-chrom CDS subset (parallel) ==="
parallel -j 24 subset_chrom "$BAKER/all.{}.vcf.gz" "$CDS/baker.{}.cds.vcf.gz" "$BED" ::: $CHROMS
echo "[$(date)] === Outgroup argentinus: per-chrom CDS subset (parallel) ==="
parallel -j 24 subset_chrom "$ARG/arg.{}.vcf.gz" "$CDS/arg.{}.cds.vcf.gz" "$BED" ::: $CHROMS

echo "[$(date)] === Concat ingroup ==="
files=(); for c in $CHROMS; do files+=("$CDS/baker.$c.cds.vcf.gz"); done
bcftools concat -Oz --threads 8 -o "$IN/baker_illecebrosus.cds.vcf.gz" "${files[@]}"
tabix -f -p vcf "$IN/baker_illecebrosus.cds.vcf.gz"

echo "[$(date)] === Concat argentinus ==="
files=(); for c in $CHROMS; do files+=("$CDS/arg.$c.cds.vcf.gz"); done
bcftools concat -Oz --threads 8 -o "$IN/argentinus.cds.vcf.gz" "${files[@]}"
tabix -f -p vcf "$IN/argentinus.cds.vcf.gz"

echo "[$(date)] === coindetii: reheader (add GT FORMAT) + index (full, already tiny) ==="
src="$ROOT/coindetti_outgroup/coindetti.snps.vcf.gz"
hdr=$(mktemp)
bcftools view -h "$src" 2>/dev/null > "$hdr"
awk '/^#CHROM/ && !d {print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; d=1} {print}' "$hdr" > "$hdr.new"
bcftools reheader -h "$hdr.new" -o "$IN/coindetii.vcf.gz" "$src"
tabix -f -p vcf "$IN/coindetii.vcf.gz"
rm -f "$hdr" "$hdr.new"

echo "[$(date)] === DONE ==="
ls -lh "$IN"/*.cds.vcf.gz
echo "variant counts:"
for f in "$IN"/baker_illecebrosus.cds.vcf.gz "$IN"/argentinus.cds.vcf.gz "$IN"/coindetii.vcf.gz; do
  echo "  $(basename $f): $(bcftools index -n $f) records"
done
