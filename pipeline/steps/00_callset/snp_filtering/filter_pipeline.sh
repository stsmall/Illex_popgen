#!/usr/bin/env bash
set -euo pipefail

# ─────────────────────────────────────────────
# INPUTS — set before running
# ─────────────────────────────────────────────
INPUT_VCF=""
INPUT_BAM=""
REFERENCE_FA=""
ACCESSIBILITY_MASK=""       # BED of accessible regions
GFF=""                      # GFF3 for snpEff database build
SNPEFF_DB=""                # snpEff database name (built from GFF, see step 3)
SNPEFF_DATA_DIR=""          # snpEff data directory
REPEATMASKER_FILE=""        # RepeatMasker .out or .gff
SAMPLE_PREFIX=""
WINDOW_SIZE=500
OUTPUT_DIR=""
THREADS=4
SCRIPTS_DIR="$(dirname "$0")"  # expects companion scripts in same directory

# ─────────────────────────────────────────────
# FILTER THRESHOLDS
# ─────────────────────────────────────────────
# Variants
VAR_MIN_QUAL=20
VAR_MIN_MQ=40
VAR_MAX_FS=60
VAR_MIN_GQ=20
VAR_MIN_DP=10
VAR_EXCHET_PARALOG=0.001
VAR_MIN_AF_PARALOG=0.1
VAR_MAX_AF_PARALOG=0.9

# Invariants (relaxed)
INV_MIN_QUAL=10
INV_MIN_MQ=30
INV_MIN_GQ=10
INV_MIN_DP=3

# ─────────────────────────────────────────────
# DERIVED PATHS
# ─────────────────────────────────────────────
mkdir -p "${OUTPUT_DIR}"/{tmp,qc,beds,annotation}

TAGGED_VCF="${OUTPUT_DIR}/tmp/00_tagged.vcf.gz"
VARIANTS_RAW="${OUTPUT_DIR}/tmp/01_variants_raw.vcf.gz"
INVARIANTS_RAW="${OUTPUT_DIR}/tmp/01_invariants_raw.vcf.gz"
VARIANTS_SNPEFF="${OUTPUT_DIR}/tmp/02_variants_snpeff.vcf.gz"
INVARIANTS_SNPEFF="${OUTPUT_DIR}/tmp/02_invariants_snpeff.vcf.gz"
HIGH_DEPTH_BED="${OUTPUT_DIR}/beds/high_depth_windows.bed"
WINDOWS_GC_BED="${OUTPUT_DIR}/beds/windows_gc.bed"
REPEAT_BED="${OUTPUT_DIR}/beds/repeats_classified.bed"
REPEAT_ANN_TAB="${OUTPUT_DIR}/annotation/repeat_annotation.tab.gz"
INACCESSIBLE_BED="${OUTPUT_DIR}/beds/inaccessible.bed"
ANN_HDR="${OUTPUT_DIR}/annotation/extra_info_header.hdr"
QC_PLOT="${OUTPUT_DIR}/qc/exchet_distribution.pdf"

# ─────────────────────────────────────────────
# STEP 1: GC-normalized depth mask
# ─────────────────────────────────────────────
echo "[1/10] Computing GC-normalized depth mask"

mosdepth \
    --by "${WINDOW_SIZE}" \
    --no-abbrev \
    --threads "${THREADS}" \
    "${OUTPUT_DIR}/tmp/${SAMPLE_PREFIX}" \
    "${INPUT_BAM}"

bedtools nuc \
    -fi "${REFERENCE_FA}" \
    -bed <(zcat "${OUTPUT_DIR}/tmp/${SAMPLE_PREFIX}.regions.bed.gz" | cut -f1-3) \
    | awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$5}' \
    > "${WINDOWS_GC_BED}"

python3 "${SCRIPTS_DIR}/gc_depth_filter.py" \
    --depth "${OUTPUT_DIR}/tmp/${SAMPLE_PREFIX}.regions.bed.gz" \
    --gc "${WINDOWS_GC_BED}" \
    --output "${HIGH_DEPTH_BED}" \
    --sd-threshold 2.0

# ─────────────────────────────────────────────
# STEP 2: Inaccessible BED
# ─────────────────────────────────────────────
echo "[2/10] Building inaccessible regions BED"

bedtools complement \
    -i "${ACCESSIBILITY_MASK}" \
    -g <(samtools fai "${REFERENCE_FA}" | cut -f1-2) \
    > "${INACCESSIBLE_BED}"

# ─────────────────────────────────────────────
# STEP 3: Build snpEff database (one-time)
# Skip if already built
# ─────────────────────────────────────────────
echo "[3/10] Building snpEff database (skip if exists)"

SNPEFF_DB_DIR="${SNPEFF_DATA_DIR}/${SNPEFF_DB}"
if [[ ! -d "${SNPEFF_DB_DIR}" ]]; then
    mkdir -p "${SNPEFF_DB_DIR}"
    cp "${GFF}" "${SNPEFF_DB_DIR}/genes.gff"
    cp "${REFERENCE_FA}" "${SNPEFF_DB_DIR}/sequences.fa"

    # Add genome entry to snpEff.config if not present
    if ! grep -q "^${SNPEFF_DB}.genome" "${SNPEFF_DATA_DIR}/../snpEff.config" 2>/dev/null; then
        echo "${SNPEFF_DB}.genome : ${SNPEFF_DB}" \
            >> "${SNPEFF_DATA_DIR}/../snpEff.config"
    fi

    snpEff build \
        -gff3 \
        -dataDir "${SNPEFF_DATA_DIR}" \
        -v "${SNPEFF_DB}"
else
    echo "    snpEff database ${SNPEFF_DB} exists, skipping build"
fi

# ─────────────────────────────────────────────
# STEP 4: Parse RepeatMasker → classified BED
#         then build bcftools-compatible annotation table
# ─────────────────────────────────────────────
echo "[4/10] Building repeat annotation table"

python3 "${SCRIPTS_DIR}/parse_repeatmasker.py" \
    --input "${REPEATMASKER_FILE}" \
    --output "${REPEAT_BED}"

# Extract all VCF positions, intersect with repeat BED
# bedtools intersect -loj keeps all VCF sites, reporting '.' if no overlap
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${INPUT_VCF}" \
| awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3, $4}' \
| bedtools intersect \
    -a stdin \
    -b "${REPEAT_BED}" \
    -loj \
| awk 'BEGIN{OFS="\t"} {
    rep = ($9 == "." ? "None" : $9);
    print $1, $3, $4, $5, rep
}' \
| bgzip > "${REPEAT_ANN_TAB}"

tabix -s1 -b2 -e2 "${REPEAT_ANN_TAB}"

# Write extra INFO header for bcftools annotate
cat > "${ANN_HDR}" << 'EOF'
##INFO=<ID=REP,Number=1,Type=String,Description="RepeatMasker class (Simple/TE/rRNA/NUMT/Other/None)">
EOF

# ─────────────────────────────────────────────
# STEP 5: Fill INFO tags
# ─────────────────────────────────────────────
echo "[5/10] Filling INFO tags"

bcftools +fill-tags "${INPUT_VCF}" \
    -Oz -o "${TAGGED_VCF}" \
    -- -t AF,HWE,ExcHet,ExpHet

bcftools index --tbi "${TAGGED_VCF}"

# ─────────────────────────────────────────────
# STEP 6: Split variants and invariants
# ─────────────────────────────────────────────
echo "[6/10] Splitting variants and invariants"

# Variants: any site with at least one alt allele observed
bcftools view -c 1 "${TAGGED_VCF}" -Oz -o "${VARIANTS_RAW}"
bcftools index --tbi "${VARIANTS_RAW}"

# Invariants: sites with no alt allele (AC=0, all samples hom-ref)
bcftools view -C 0 "${TAGGED_VCF}" -Oz -o "${INVARIANTS_RAW}"
bcftools index --tbi "${INVARIANTS_RAW}"

# ─────────────────────────────────────────────
# STEP 7: QC diagnostic on variants
# ─────────────────────────────────────────────
echo "[7/10] Generating QC plots"

bcftools query \
    -i 'AF > 0.1 && AF < 0.9' \
    -f '%CHROM\t%POS\t%AF\t%ExcHet\t%INFO/DP\n' \
    "${VARIANTS_RAW}" \
| Rscript - "${QC_PLOT}" << 'EOF'
args <- commandArgs(trailingOnly=TRUE)
outpdf <- args[1]
d <- read.table("stdin", header=FALSE,
                col.names=c("chr","pos","af","exchet","dp"))
pdf(outpdf)
par(mfrow=c(1,2))
hist(-log10(d$exchet), breaks=100,
     main="ExcHet (mid-AF sites)", xlab="-log10(ExcHet)")
abline(v=3, col="red", lty=2)
legend("topright", legend="p=0.001", col="red", lty=2)
plot(d$af, -log10(d$exchet), pch=20, cex=0.3,
     xlab="AF", ylab="-log10(ExcHet)", main="ExcHet vs AF")
abline(h=3, col="red", lty=2)
dev.off()
cat("QC plot written to:", outpdf, "\n")
EOF

echo "    >>> Review ${QC_PLOT} and adjust VAR_EXCHET_PARALOG if needed <<<"

# ─────────────────────────────────────────────
# STEP 8: snpEff annotation
# Run on variants and invariants separately
# ─────────────────────────────────────────────
echo "[8/10] Running snpEff annotation"

snpEff ann \
    -dataDir "${SNPEFF_DATA_DIR}" \
    -noStats \
    -noLog \
    "${SNPEFF_DB}" \
    "${VARIANTS_RAW}" \
| bgzip > "${VARIANTS_SNPEFF}"
bcftools index --tbi "${VARIANTS_SNPEFF}"

snpEff ann \
    -dataDir "${SNPEFF_DATA_DIR}" \
    -noStats \
    -noLog \
    "${SNPEFF_DB}" \
    "${INVARIANTS_RAW}" \
| bgzip > "${INVARIANTS_SNPEFF}"
bcftools index --tbi "${INVARIANTS_SNPEFF}"

# ─────────────────────────────────────────────
# STEP 9: Add repeat annotation + filter + produce outputs
# ─────────────────────────────────────────────
echo "[9/10] Adding repeat annotation and filtering"

filter_and_annotate() {
    local INPUT="$1"
    local MODE="$2"          # "variant" or "invariant"
    local SOFTFILTER_OUT="$3"
    local HARDFILT_OUT="$4"

    local REP_ANNOTATED="${OUTPUT_DIR}/tmp/${MODE}_rep_annotated.vcf.gz"

    bcftools annotate \
        -a "${REPEAT_ANN_TAB}" \
        -c CHROM,POS,REF,ALT,INFO/REP \
        -h "${ANN_HDR}" \
        "${INPUT}" \
        -Oz -o "${REP_ANNOTATED}"
    bcftools index --tbi "${REP_ANNOTATED}"

    if [[ "${MODE}" == "variant" ]]; then
        bcftools norm -m -any "${REP_ANNOTATED}" \
        | bcftools filter \
            --soft-filter LowSiteQual \
            -m + \
            -e "QUAL < ${VAR_MIN_QUAL} || MQ < ${VAR_MIN_MQ} || FS > ${VAR_MAX_FS}" \
        | bcftools filter \
            --soft-filter PARALOG \
            -m + \
            -e "AF > ${VAR_MIN_AF_PARALOG} && AF < ${VAR_MAX_AF_PARALOG} && ExcHet < ${VAR_EXCHET_PARALOG}" \
        | bcftools filter \
            --soft-filter LowGQ \
            -m + \
            -e "FMT/GQ < ${VAR_MIN_GQ}" \
        | bcftools filter \
            --soft-filter LowDP \
            -m + \
            -e "FMT/DP < ${VAR_MIN_DP}" \
        | bcftools filter \
            --soft-filter BadAB \
            -m + \
            -e '(GT="het") & ((FMT/AD[*:1]/FMT/DP < 0.3) | (FMT/AD[*:1]/FMT/DP > 0.7))' \
        | bcftools filter \
            --soft-filter HIGH_DEPTH_GC \
            -m + \
            -T "${HIGH_DEPTH_BED}" \
        | bcftools filter \
            --soft-filter INACCESSIBLE \
            -m + \
            -T "${INACCESSIBLE_BED}" \
            -Oz -o "${SOFTFILTER_OUT}"
    else
        # Invariant: relaxed thresholds, no AB/FS/paralog filters
        bcftools filter \
            --soft-filter LowSiteQual \
            -m + \
            -e "QUAL < ${INV_MIN_QUAL} || MQ < ${INV_MIN_MQ}" \
            "${REP_ANNOTATED}" \
        | bcftools filter \
            --soft-filter LowGQ \
            -m + \
            -e "FMT/GQ < ${INV_MIN_GQ}" \
        | bcftools filter \
            --soft-filter LowDP \
            -m + \
            -e "FMT/DP < ${INV_MIN_DP}" \
        | bcftools filter \
            --soft-filter HIGH_DEPTH_GC \
            -m + \
            -T "${HIGH_DEPTH_BED}" \
        | bcftools filter \
            --soft-filter INACCESSIBLE \
            -m + \
            -T "${INACCESSIBLE_BED}" \
            -Oz -o "${SOFTFILTER_OUT}"
    fi

    bcftools index --tbi "${SOFTFILTER_OUT}"

    # Hard filter: PASS only
    bcftools view -f PASS \
        "${SOFTFILTER_OUT}" \
        -Oz -o "${HARDFILT_OUT}"
    bcftools index --tbi "${HARDFILT_OUT}"

    # Companion TSV via SnpSift for clean ANN field parsing
    SnpSift extractFields \
        "${HARDFILT_OUT}" \
        CHROM POS REF ALT \
        "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].BIOTYPE" \
        INFO/REP \
        > "${HARDFILT_OUT%.vcf.gz}.tsv"
}

# ─────────────────────────────────────────────
# STEP 10: Run filter_and_annotate for both streams
# ─────────────────────────────────────────────
echo "[10/10] Producing final outputs"

filter_and_annotate \
    "${VARIANTS_SNPEFF}" \
    "variant" \
    "${OUTPUT_DIR}/variants_softfilter.vcf.gz" \
    "${OUTPUT_DIR}/variants_filt.vcf.gz"

filter_and_annotate \
    "${INVARIANTS_SNPEFF}" \
    "invariant" \
    "${OUTPUT_DIR}/invariant_softfilter.vcf.gz" \
    "${OUTPUT_DIR}/invariant_filt.vcf.gz"

# ─────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────
echo "Done. Filter summary:"

{
    for VCF in variants_filt variants_softfilter invariant_filt invariant_softfilter; do
        echo "--- ${VCF} ---"
        bcftools stats "${OUTPUT_DIR}/${VCF}.vcf.gz" | grep "^SN"
        bcftools query -f '%FILTER\n' "${OUTPUT_DIR}/${VCF}.vcf.gz" \
            | sort | uniq -c | sort -rn
        echo ""
    done
} > "${OUTPUT_DIR}/qc/filter_summary.txt"

cat "${OUTPUT_DIR}/qc/filter_summary.txt"
echo "All outputs written to ${OUTPUT_DIR}"
