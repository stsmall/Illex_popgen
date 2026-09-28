#!/usr/bin/env bash
# build_shared_inputs.sh
# Compute chrom-independent shared filter inputs ONCE (genome-wide).
# Outputs to: steps/00_callset/filtered/_shared/
#
# Shared outputs produced here:
#   _shared/beds/high_depth_windows.bed   — GC-normalized high-depth mask (mosdepth)
#   _shared/beds/inaccessible.bed         — complement of $ACCESS vs genome .fai
#   _shared/annotation/extra_info_header.hdr — bcftools annotate header for REP tag
#
# REP annotation is handled per-chrom in run_filter.sh via bedtools intersect
# against repeats.merged.bed at runtime (no pre-built genome-wide table needed).
# This avoids processing ~38M chr-Z positions x 46 chroms upfront.
#
# NOTE on GC-depth proxy: mosdepth is run on a SINGLE representative BAM
# (SQ233K001.bam). This is a single-sample proxy; the mask reflects per-sample
# depth variation rather than a multi-sample mean. Limitation: may miss
# cross-sample high-depth regions; may include sample-specific artefacts.
# Documented per task brief.

set -euo pipefail

# ─────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh

SHARED_DIR="$ANA/steps/00_callset/filtered/_shared"
SCRIPTS_DIR="$ANA/steps/00_callset/snp_filtering"
INPUT_BAM=/sietch_colab/data_share/illex/popgen_data/seq_data/baker_2025/mapping/dedup/SQ233K001.bam
WINDOW_SIZE=500
THREADS=16

# mosdepth and gc_depth_filter.py use bioinfo-buddy env (has mosdepth + python3 w/ pandas/numpy/statsmodels)
# bcftools/bedtools/bgzip/tabix use ENV_BCF (mkado-vcf) for VCF operations
ENV_MOSDEEP=/home/ssmall/miniforge3/envs/bioinfo-buddy

# ─────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────
mkdir -p "$SHARED_DIR"/{tmp,beds,annotation}

HIGH_DEPTH_BED="$SHARED_DIR/beds/high_depth_windows.bed"
WINDOWS_GC_BED="$SHARED_DIR/beds/windows_gc.bed"
INACCESSIBLE_BED="$SHARED_DIR/beds/inaccessible.bed"
ANN_HDR="$SHARED_DIR/annotation/extra_info_header.hdr"
MOSDEPTH_PREFIX="$SHARED_DIR/tmp/mosdepth_proxy"

# ─────────────────────────────────────────────
# STEP A: GC-normalized high-depth mask
# ─────────────────────────────────────────────
echo "[A/B] Computing GC-normalized depth mask (genome-wide, single-sample proxy)"
echo "      BAM: ${INPUT_BAM}"

if [[ ! -f "${MOSDEPTH_PREFIX}.regions.bed.gz" ]]; then
    mamba run --no-banner -p "$ENV_MOSDEEP" \
        mosdepth \
            --by "${WINDOW_SIZE}" \
            --threads "${THREADS}" \
            "${MOSDEPTH_PREFIX}" \
            "${INPUT_BAM}"
    echo "  mosdepth done"
else
    echo "  mosdepth output exists, skipping"
fi

if [[ ! -f "${WINDOWS_GC_BED}" ]]; then
    mamba run --no-banner -p "$ENV_MOSDEEP" bash -c "
        bedtools nuc \
            -fi '${REF}' \
            -bed <(zcat '${MOSDEPTH_PREFIX}.regions.bed.gz' | cut -f1-3) \
        | awk 'NR>1 {print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$5}' \
        > '${WINDOWS_GC_BED}'
    "
    echo "  GC content computed: ${WINDOWS_GC_BED}"
else
    echo "  GC content BED exists, skipping"
fi

if [[ ! -f "${HIGH_DEPTH_BED}" ]]; then
    mamba run --no-banner -p "$ENV_MOSDEEP" \
        python3 "${SCRIPTS_DIR}/gc_depth_filter.py" \
            --depth "${MOSDEPTH_PREFIX}.regions.bed.gz" \
            --gc "${WINDOWS_GC_BED}" \
            --output "${HIGH_DEPTH_BED}" \
            --sd-threshold 2.0
    echo "  High-depth BED written: ${HIGH_DEPTH_BED}"
else
    echo "  High-depth BED exists, skipping"
fi

# ─────────────────────────────────────────────
# STEP B: Inaccessible BED (complement of ACCESS)
# ─────────────────────────────────────────────
echo "[B/B] Building inaccessible regions BED"

if [[ ! -f "${INACCESSIBLE_BED}" ]]; then
    # Sort ACCESS BED to match FAI chromosome order (FAI is 1..45,Z; ACCESS may differ)
    # Use sortBed with -faidx to enforce FAI order
    SORTED_ACCESS="$SHARED_DIR/tmp/accessible_sorted.bed"
    mamba run --no-banner -p "$ENV_MOSDEEP" bash -c "
        bedtools sort -i '${ACCESS}' -faidx '${REF}.fai' > '${SORTED_ACCESS}'
        bedtools complement \
            -i '${SORTED_ACCESS}' \
            -g <(cut -f1-2 '${REF}.fai') \
            > '${INACCESSIBLE_BED}'
    "
    echo "  Inaccessible BED written: ${INACCESSIBLE_BED}"
else
    echo "  Inaccessible BED exists, skipping"
fi

# ─────────────────────────────────────────────
# Write REP annotation header (used by run_filter.sh per-chrom)
# ─────────────────────────────────────────────
cat > "${ANN_HDR}" << 'EOF'
##INFO=<ID=REP,Number=1,Type=String,Description="Repeat overlap (from repeats.merged.bed, unclassified merged union of EarlGrey+TRASH): Repeat or None">
EOF
echo "  REP header written: ${ANN_HDR}"

echo ""
echo "=== Shared inputs complete ==="
echo "  HIGH_DEPTH_BED:   ${HIGH_DEPTH_BED}"
echo "  INACCESSIBLE_BED: ${INACCESSIBLE_BED}"
echo "  ANN_HDR:          ${ANN_HDR}"
echo ""
echo "NOTE: REP annotation table is built per-chrom in run_filter.sh"
echo "      (bedtools intersect against repeats.merged.bed at runtime)."
