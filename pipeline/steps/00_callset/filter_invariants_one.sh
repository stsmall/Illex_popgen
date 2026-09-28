#!/usr/bin/env bash
# filter_invariants_one.sh — per-chrom INVARIANT filter driver (Task 0.2 re-arch)
#
# Usage:
#   CHR=Z bash filter_invariants_one.sh
#   CHR=Z REGION=Z:1-3000000 bash filter_invariants_one.sh   # slice mode for testing
#
# Outputs (filtered/<chr>/):
#   invariant_softfilter.vcf.gz  — all invariants, FILTER column populated
#   invariant_filt.vcf.gz        — PASS-only
#
# Design: MINIMAL and STREAMED
#   - NO snpEff (annotation not meaningful for invariant sites)
#   - fill-tags computes ONLY F_MISSING (cheap; for the site call-rate filter)
#   - NO REP annotation (not needed for invariant callable sites)
#   - GT-only callset: no FMT/GQ,FMT/DP -> those filters removed
#   - One read of the merged VCF, streamed through bcftools filter chain
#   - Only TWO writes: softfilter_out (with FILTER column) + filt_out (PASS-only)
#   - Relaxed thresholds vs variants (per task brief defaults)
#
# Key envs (per interface-notes §7):
#   bioinfo-buddy  — bcftools 1.22 + bgzip/tabix
#
# BUG FIX vs run_filter.sh:
#   Original used -T <bed> in bcftools filter, which is --targets-file and SUBSETS
#   output to BED sites only (drops everything outside the BED).
#   Correct flag for soft-filtering regions: --mask-file <bed>
#   Verified: -T test: 101 records in, 54 out (subset); --mask-file: 101 out (all preserved).

set -euo pipefail

# ─── Config ────────────────────────────────────────────────────────────────────
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh

[[ -z "${CHR:-}" ]] && { echo "ERROR: CHR must be set (e.g. CHR=Z bash $0)" >&2; exit 1; }

SHARED_DIR="$ANA/steps/00_callset/filtered/_shared"
OUTPUT_DIR="$ANA/steps/00_callset/filtered/${CHR}"
INPUT_VCF="$ANA/steps/00_callset/merged/all.${CHR}.vcf.gz"

ENV_BUDDY=/home/ssmall/miniforge3/envs/bioinfo-buddy
THREADS=16

# Filter thresholds — invariants (GT-only callset: no FMT/GQ,FMT/DP available.
# Relaxed site QUAL/MQ + site call-rate F_MISSING + region masks.)
INV_MIN_QUAL=10
INV_MIN_MQ=30
INV_MAX_FMISSING=0.5

HIGH_DEPTH_BED="$SHARED_DIR/beds/high_depth_windows.bed"
INACCESSIBLE_BED="$SHARED_DIR/beds/inaccessible.bed"

# Optional region for slice testing (e.g. REGION=Z:1-3000000)
REGION="${REGION:-}"

# ─── Validate inputs ────────────────────────────────────────────────────────────
for f in "$HIGH_DEPTH_BED" "$INACCESSIBLE_BED" "$INPUT_VCF"; do
    [[ -f "$f" ]] || { echo "ERROR: missing required file: $f" >&2; exit 1; }
done

# ─── Setup dirs ─────────────────────────────────────────────────────────────────
mkdir -p "$OUTPUT_DIR"/qc

echo "============================================================"
echo "filter_invariants_one.sh  CHR=${CHR}${REGION:+  REGION=${REGION}}"
echo "INPUT:  ${INPUT_VCF}"
echo "OUTPUT: ${OUTPUT_DIR}"
echo "============================================================"

# ─── Single-pass streamed filter chain ─────────────────────────────────────────
# bcftools view -C0 → streamed invariant extraction
# chained bcftools filter with -Ou between stages (uncompressed pipe, no intermediate writes)
# -Oz only at the final write
#
# CRITICAL: --mask-file (not -T) for region tagging — see header note.
#
# Two writes total:
#   1. invariant_softfilter.vcf.gz  (all invariants with FILTER column)
#   2. invariant_filt.vcf.gz        (PASS-only)

SOFTFILTER_OUT="$OUTPUT_DIR/invariant_softfilter.vcf.gz"
HARDFILT_OUT="$OUTPUT_DIR/invariant_filt.vcf.gz"

# DATA-INTEGRITY: same corrupt ALT=* records as the variant driver crash bcftools view -C0
# (see report §0). An inline tabix|awk sanitizer is too slow at genome scale. RESOLUTION:
# fix ONCE upstream; if all.<chr>.clean.vcf.gz exists it is preferred (fast indexed path).
CLEAN_VCF="$ANA/steps/00_callset/merged/all.${CHR}.clean.vcf.gz"
if [[ -f "$CLEAN_VCF" ]]; then
    SRC_VCF="$CLEAN_VCF"
    echo "  Using sanitized input: $CLEAN_VCF"
else
    SRC_VCF="$INPUT_VCF"
fi

REGION_ARG=""
[[ -n "$REGION" ]] && REGION_ARG="-r $REGION"

echo "[1/2] Streaming invariant filter chain..."
if [[ ! -f "$SOFTFILTER_OUT" ]]; then
    T_START=$(date +%s)

    mamba run --no-banner -p "$ENV_BUDDY" bash -c "
        bcftools view -C0 ${REGION_ARG} -Ou '${SRC_VCF}' \
        | bcftools +fill-tags -Ou -- -t F_MISSING \
        | bcftools filter \
            --soft-filter LowSiteQual \
            -m+ \
            -e 'QUAL < ${INV_MIN_QUAL} || MQ < ${INV_MIN_MQ}' \
            -Ou \
        | bcftools filter \
            --soft-filter MISSING \
            -m+ \
            -e 'F_MISSING >= ${INV_MAX_FMISSING}' \
            -Ou \
        | bcftools filter \
            --soft-filter HIGH_DEPTH_GC \
            -m+ \
            --mask-file '${HIGH_DEPTH_BED}' \
            -Ou \
        | bcftools filter \
            --soft-filter INACCESSIBLE \
            -m+ \
            --mask-file '${INACCESSIBLE_BED}' \
            -Oz -o '${SOFTFILTER_OUT}'
    "
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$SOFTFILTER_OUT"

    T_END=$(date +%s)
    T_ELAPSED=$(( T_END - T_START ))

    # Report throughput for the slice (for full-genome concurrency estimation)
    N_RECORDS=$(mamba run --no-banner -p "$ENV_BUDDY" bcftools stats "$SOFTFILTER_OUT" | awk -F'\t' '/^SN.*records/{print $4}')
    if [[ -n "$REGION" && "$T_ELAPSED" -gt 0 && -n "${N_RECORDS:-}" ]]; then
        SITES_PER_MIN=$(( N_RECORDS * 60 / T_ELAPSED ))
        echo "  Throughput: ${N_RECORDS} records in ${T_ELAPSED}s = ~${SITES_PER_MIN} sites/min"
    fi
    echo "  Softfilter VCF written: ${N_RECORDS:-unknown} records"
else
    echo "  Softfilter VCF exists, skipping."
fi

echo "[2/2] Extracting PASS-only invariants..."
if [[ ! -f "$HARDFILT_OUT" ]]; then
    mamba run --no-banner -p "$ENV_BUDDY" \
        bcftools view -f PASS "$SOFTFILTER_OUT" -Oz -o "$HARDFILT_OUT"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$HARDFILT_OUT"
    PASS_N=$(mamba run --no-banner -p "$ENV_BUDDY" bcftools stats "$HARDFILT_OUT" | awk -F'\t' '/^SN.*records/{print $4}')
    echo "  PASS-only VCF: ${PASS_N} records"
else
    echo "  PASS VCF exists, skipping."
fi

# ─── Summary ────────────────────────────────────────────────────────────────────
echo ""
echo "=== Invariant filter summary CHR=${CHR} $(date) ==="
{
    echo "=== CHR=${CHR} invariant filter summary $(date) ==="
    for label in invariant_filt invariant_softfilter; do
        VCF="$OUTPUT_DIR/${label}.vcf.gz"
        [[ -f "$VCF" ]] || continue
        echo "--- ${label} ---"
        mamba run --no-banner -p "$ENV_BUDDY" bcftools stats "$VCF" | grep "^SN"
        mamba run --no-banner -p "$ENV_BUDDY" bcftools query -f '%FILTER\n' "$VCF" \
            | sort | uniq -c | sort -rn
        echo ""
    done
} | tee "$OUTPUT_DIR/qc/invariant_filter_summary.txt"

echo "Done. Outputs in: ${OUTPUT_DIR}"
