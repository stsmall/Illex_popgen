#!/usr/bin/env bash
# run_filter.sh — per-chrom filtering driver for snp_filtering pipeline
#
# Usage:
#   CHR=Z bash run_filter.sh
#   CHR=1 bash run_filter.sh
#
# Requires shared inputs to be built first:
#   bash build_shared_inputs.sh
#
# Shared inputs consumed:
#   filtered/_shared/beds/high_depth_windows.bed
#   filtered/_shared/beds/inaccessible.bed
#   filtered/_shared/annotation/extra_info_header.hdr
#
# REP annotation strategy (no usable RepeatMasker .out/.gff for primary assembly):
#   Annotate REP INFO field via bedtools intersect against repeats.merged.bed
#   (3-col BED, merged union of EarlGrey + TRASH outputs, chroms 1-45 and Z).
#   All overlaps classified as REP=Repeat; non-overlapping sites get REP=None.
#   parse_repeatmasker.py is NOT used (no RepeatMasker .out input available).
#
# snpEff DB: pre-built as 'illex' (build step skipped per interface-notes §5).
# Java >=23 is required for snpEff 5.4c; mkado-vcf env has Java 11 only —
# we use the java from mkado_env which has Java 25.
#
# Filter thresholds (exact per brief/README defaults):
#   Variants:   QUAL>=20, MQ>=40, FS<=60, GQ>=20/sample, DP>=10/sample,
#               ExcHet<0.001 at 0.1<AF<0.9 (PARALOG), BadAB, HIGH_DEPTH_GC, INACCESSIBLE
#   Invariants: QUAL>=10, MQ>=30, GQ>=10/sample, DP>=3/sample,
#               HIGH_DEPTH_GC, INACCESSIBLE

set -euo pipefail

# ─────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh

if [[ -z "${CHR:-}" ]]; then
    echo "ERROR: CHR must be set (e.g. CHR=Z bash run_filter.sh)" >&2
    exit 1
fi

SHARED_DIR="$ANA/steps/00_callset/filtered/_shared"
OUTPUT_DIR="$ANA/steps/00_callset/filtered/${CHR}"
INPUT_VCF="$ANA/steps/00_callset/merged/all.${CHR}.vcf.gz"
REP_BED=/sietch_colab/data_share/illex/annotations/repeat_masking_dir/temp_/repeats.merged.bed
SCRIPTS_DIR="$ANA/steps/00_callset/snp_filtering"

# snpEff settings (DB pre-built; build step skipped)
SNPEFF_DB="illex"
SNPEFF_DATA_DIR="/home/ssmall/programs/snpEff/data"
SNPEFF_JAR="/home/ssmall/programs/snpEff/snpEff.jar"
SNPSIFT_JAR="/home/ssmall/programs/snpEff/SnpSift.jar"
JAVA23="/home/ssmall/miniforge3/envs/mkado_env/lib/jvm/bin/java"

THREADS=16

# bioinfo-buddy has bcftools 1.22 with fill-tags plugin + bedtools + bgzip/tabix + python3 w/ deps
# mkado-vcf uses /home/ssmall/bin/bcftools 1.21 which lacks plugin dir — don't use for VCF ops
ENV_BUDDY=/home/ssmall/miniforge3/envs/bioinfo-buddy

# Filter thresholds — variants
VAR_MIN_QUAL=20
VAR_MIN_MQ=40
VAR_MAX_FS=60
VAR_MIN_GQ=20
VAR_MIN_DP=10
VAR_EXCHET_PARALOG=0.001
VAR_MIN_AF_PARALOG=0.1
VAR_MAX_AF_PARALOG=0.9

# Filter thresholds — invariants (relaxed)
INV_MIN_QUAL=10
INV_MIN_MQ=30
INV_MIN_GQ=10
INV_MIN_DP=3

# Shared inputs
HIGH_DEPTH_BED="$SHARED_DIR/beds/high_depth_windows.bed"
INACCESSIBLE_BED="$SHARED_DIR/beds/inaccessible.bed"
ANN_HDR="$SHARED_DIR/annotation/extra_info_header.hdr"

# ─────────────────────────────────────────────
# Validate shared inputs exist
# ─────────────────────────────────────────────
for f in "$HIGH_DEPTH_BED" "$INACCESSIBLE_BED" "$ANN_HDR" "$INPUT_VCF"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file missing: $f" >&2
        echo "Run build_shared_inputs.sh first." >&2
        exit 1
    fi
done

# ─────────────────────────────────────────────
# Setup output directories
# ─────────────────────────────────────────────
mkdir -p "$OUTPUT_DIR"/{tmp,qc,beds,annotation}

# ─────────────────────────────────────────────
# Concurrency guard: refuse to start if another run_filter.sh is
# already processing this same chrom (prevents two processes writing
# the same intermediate files and corrupting them).
# ─────────────────────────────────────────────
LOCKFILE="$OUTPUT_DIR/.run_filter.lock"
if [[ -e "$LOCKFILE" ]]; then
    OTHER_PID=$(cat "$LOCKFILE" 2>/dev/null || echo "")
    if [[ -n "$OTHER_PID" ]] && kill -0 "$OTHER_PID" 2>/dev/null; then
        echo "ERROR: run_filter.sh already running for CHR=${CHR} (PID $OTHER_PID). Refusing to start a second run." >&2
        exit 1
    fi
    echo "  Stale lockfile found (PID ${OTHER_PID:-?} not running) — reclaiming."
fi
echo "$$" > "$LOCKFILE"
trap 'rm -f "$LOCKFILE"' EXIT

SAMPLE_PREFIX="illex_${CHR}"

# Intermediate files
TAGGED_VCF="$OUTPUT_DIR/tmp/00_tagged.vcf.gz"
VARIANTS_RAW="$OUTPUT_DIR/tmp/01_variants_raw.vcf.gz"
INVARIANTS_RAW="$OUTPUT_DIR/tmp/01_invariants_raw.vcf.gz"
VARIANTS_SNPEFF="$OUTPUT_DIR/tmp/02_variants_snpeff.vcf.gz"
INVARIANTS_SNPEFF="$OUTPUT_DIR/tmp/02_invariants_snpeff.vcf.gz"
REPEAT_ANN_TAB="$OUTPUT_DIR/annotation/repeat_annotation.tab.gz"
QC_PLOT="$OUTPUT_DIR/qc/exchet_distribution.pdf"

echo "============================================================"
echo "run_filter.sh  CHR=${CHR}"
echo "INPUT:  ${INPUT_VCF}"
echo "OUTPUT: ${OUTPUT_DIR}"
echo "============================================================"

# ─────────────────────────────────────────────
# STEP 5: Fill INFO tags (AF, HWE, ExcHet, ExpHet)
# ─────────────────────────────────────────────
echo "[5/10] Filling INFO tags (AF, HWE, ExcHet, ExpHet)"

if [[ ! -f "$TAGGED_VCF" ]]; then
    # ExpHet not supported in bcftools 1.22 fill-tags; dropped (not used in any filter)
    mamba run --no-banner -p "$ENV_BUDDY" \
        bcftools +fill-tags "$INPUT_VCF" \
            -Oz -o "$TAGGED_VCF" \
            -- -t AF,HWE,ExcHet
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$TAGGED_VCF"
    echo "  Tagged VCF written."
else
    echo "  Tagged VCF exists, skipping."
fi

# ─────────────────────────────────────────────
# STEP 6: Split variants and invariants
# ─────────────────────────────────────────────
echo "[6/10] Splitting variants and invariants"

if [[ ! -f "$VARIANTS_RAW" ]]; then
    # Variants: at least one alt allele observed (-c 1)
    mamba run --no-banner -p "$ENV_BUDDY" \
        bcftools view -c 1 "$TAGGED_VCF" -Oz -o "$VARIANTS_RAW"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$VARIANTS_RAW"
    echo "  Variants raw written."
else
    echo "  Variants raw exists, skipping."
fi

if [[ ! -f "$INVARIANTS_RAW" ]]; then
    # Invariants: no alt allele (AC=0, all hom-ref) (-C 0)
    mamba run --no-banner -p "$ENV_BUDDY" \
        bcftools view -C 0 "$TAGGED_VCF" -Oz -o "$INVARIANTS_RAW"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$INVARIANTS_RAW"
    echo "  Invariants raw written."
else
    echo "  Invariants raw exists, skipping."
fi

# ─────────────────────────────────────────────
# STEP 7: QC diagnostic — ExcHet distribution
# ─────────────────────────────────────────────
echo "[7/10] Generating ExcHet QC plot"

if [[ ! -f "$QC_PLOT" ]]; then
    # Write R script to a temp file (can't combine pipe + heredoc for Rscript stdin)
    QC_RSCRIPT="$OUTPUT_DIR/tmp/exchet_qc.R"
    cat > "$QC_RSCRIPT" << 'REOF'
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outpdf <- args[2]
d <- read.table(infile, header=FALSE,
                col.names=c("chr","pos","af","exchet","dp"))
cat("Mid-AF sites (0.1<AF<0.9):", nrow(d), "\n")
if (nrow(d) == 0) { stop("No mid-AF sites found for QC plot.") }

# Stats
n_below_thresh <- sum(d$exchet < 0.001, na.rm=TRUE)
cat("Sites with ExcHet < 0.001:", n_below_thresh, "(", round(100*n_below_thresh/nrow(d),2), "%)\n")
cat("ExcHet quantiles (mid-AF sites):\n")
print(quantile(d$exchet, c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE))

pdf(outpdf)
par(mfrow=c(1,2))
hist(-log10(d$exchet), breaks=100,
     main=paste0("ExcHet (mid-AF sites, n=", nrow(d), ")"),
     xlab="-log10(ExcHet)")
abline(v=3, col="red", lty=2)
legend("topright", legend="p=0.001", col="red", lty=2)
plot(d$af, -log10(d$exchet), pch=20, cex=0.3,
     xlab="AF", ylab="-log10(ExcHet)", main="ExcHet vs AF")
abline(h=3, col="red", lty=2)
dev.off()
cat("QC plot written to:", outpdf, "\n")
REOF

    QC_DATA="$OUTPUT_DIR/tmp/exchet_data.tsv"
    mamba run --no-banner -p "$ENV_BUDDY" \
        bcftools query \
            -i 'AF > 0.1 && AF < 0.9' \
            -f '%CHROM\t%POS\t%AF\t%ExcHet\t%INFO/DP\n' \
            "$VARIANTS_RAW" \
        > "$QC_DATA"

    Rscript "$QC_RSCRIPT" "$QC_DATA" "$QC_PLOT"
    echo "  QC plot: $QC_PLOT"
else
    echo "  QC plot exists, skipping."
fi

# ─────────────────────────────────────────────
# STEP 8: snpEff annotation (DB pre-built; build skipped)
# ─────────────────────────────────────────────
echo "[8/10] Running snpEff annotation (DB: ${SNPEFF_DB})"

if [[ ! -f "$VARIANTS_SNPEFF" ]]; then
    "$JAVA23" -Xmx16g -jar "$SNPEFF_JAR" ann \
        -dataDir "$SNPEFF_DATA_DIR" \
        -noStats \
        -noLog \
        -noCheckCds \
        -noCheckProtein \
        "$SNPEFF_DB" \
        "$VARIANTS_RAW" \
    | mamba run --no-banner -p "$ENV_BUDDY" bgzip > "$VARIANTS_SNPEFF"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$VARIANTS_SNPEFF"
    echo "  Variants snpEff done."
else
    echo "  Variants snpEff exists, skipping."
fi

if [[ ! -f "$INVARIANTS_SNPEFF" ]]; then
    "$JAVA23" -Xmx16g -jar "$SNPEFF_JAR" ann \
        -dataDir "$SNPEFF_DATA_DIR" \
        -noStats \
        -noLog \
        -noCheckCds \
        -noCheckProtein \
        "$SNPEFF_DB" \
        "$INVARIANTS_RAW" \
    | mamba run --no-banner -p "$ENV_BUDDY" bgzip > "$INVARIANTS_SNPEFF"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$INVARIANTS_SNPEFF"
    echo "  Invariants snpEff done."
else
    echo "  Invariants snpEff exists, skipping."
fi

# ─────────────────────────────────────────────
# STEP 4 (done per-chrom): Build REP annotation table for this chrom
# Strategy: bedtools intersect VCF positions against repeats.merged.bed
# Output: bgzipped CHROM/POS/REF/ALT/REP table for bcftools annotate
# All overlaps = REP=Repeat; non-overlapping = REP=None (binary flag)
# ─────────────────────────────────────────────
echo "[4/10] Building REP annotation table for chr ${CHR}"

if [[ ! -f "$REPEAT_ANN_TAB" ]]; then
    # Extract all positions from the chr-specific tagged VCF
    # Build BED5, intersect -loj with repeat BED, output annotation table
    mamba run --no-banner -p "$ENV_BUDDY" bash -c "
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' '${TAGGED_VCF}' \
        | awk 'BEGIN{OFS=\"\t\"} {print \$1, \$2-1, \$2, \$3, \$4}' \
        | bedtools intersect \
            -a stdin \
            -b '${REP_BED}' \
            -loj \
        | awk 'BEGIN{OFS=\"\t\"} {
            rep = (\$6 == \".\" ? \"None\" : \"Repeat\");
            print \$1, \$3, \$4, \$5, rep
        }' \
        | bgzip > '${REPEAT_ANN_TAB}'
        tabix -s1 -b2 -e2 '${REPEAT_ANN_TAB}'
    "
    echo "  REP annotation table built: ${REPEAT_ANN_TAB}"
else
    echo "  REP annotation table exists, skipping."
fi

# ─────────────────────────────────────────────
# STEP 9 + 10: Add REP annotation, apply filters, produce outputs
# ─────────────────────────────────────────────
echo "[9-10/10] Adding REP annotation and filtering"

filter_and_annotate() {
    local INPUT="$1"
    local MODE="$2"           # "variant" or "invariant"
    local SOFTFILTER_OUT="$3"
    local HARDFILT_OUT="$4"

    local REP_ANNOTATED="$OUTPUT_DIR/tmp/${MODE}_rep_annotated.vcf.gz"

    echo "  [${MODE}] Adding REP annotation..."
    if [[ ! -f "$REP_ANNOTATED" ]]; then
        mamba run --no-banner -p "$ENV_BUDDY" \
            bcftools annotate \
                -a "$REPEAT_ANN_TAB" \
                -c CHROM,POS,REF,ALT,INFO/REP \
                -h "$ANN_HDR" \
                "$INPUT" \
                -Oz -o "$REP_ANNOTATED"
        mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$REP_ANNOTATED"
    else
        echo "  [${MODE}] REP-annotated exists, skipping annotation."
    fi

    echo "  [${MODE}] Applying filters → softfilter..."
    if [[ ! -f "$SOFTFILTER_OUT" ]]; then
        if [[ "${MODE}" == "variant" ]]; then
            mamba run --no-banner -p "$ENV_BUDDY" bash -c "
                bcftools norm -m -any '${REP_ANNOTATED}' \
                | bcftools filter \
                    --soft-filter LowSiteQual \
                    -m + \
                    -e 'QUAL < ${VAR_MIN_QUAL} || MQ < ${VAR_MIN_MQ} || FS > ${VAR_MAX_FS}' \
                | bcftools filter \
                    --soft-filter PARALOG \
                    -m + \
                    -e 'AF > ${VAR_MIN_AF_PARALOG} && AF < ${VAR_MAX_AF_PARALOG} && ExcHet < ${VAR_EXCHET_PARALOG}' \
                | bcftools filter \
                    --soft-filter LowGQ \
                    -m + \
                    -e 'FMT/GQ < ${VAR_MIN_GQ}' \
                | bcftools filter \
                    --soft-filter LowDP \
                    -m + \
                    -e 'FMT/DP < ${VAR_MIN_DP}' \
                | bcftools filter \
                    --soft-filter BadAB \
                    -m + \
                    -e '(GT=\"het\") & ((FMT/AD[*:1]/FMT/DP < 0.3) | (FMT/AD[*:1]/FMT/DP > 0.7))' \
                | bcftools filter \
                    --soft-filter HIGH_DEPTH_GC \
                    -m + \
                    -T '${HIGH_DEPTH_BED}' \
                | bcftools filter \
                    --soft-filter INACCESSIBLE \
                    -m + \
                    -T '${INACCESSIBLE_BED}' \
                    -Oz -o '${SOFTFILTER_OUT}'
            "
        else
            # Invariant: relaxed thresholds, no AB/FS/paralog filters
            mamba run --no-banner -p "$ENV_BUDDY" bash -c "
                bcftools filter \
                    --soft-filter LowSiteQual \
                    -m + \
                    -e 'QUAL < ${INV_MIN_QUAL} || MQ < ${INV_MIN_MQ}' \
                    '${REP_ANNOTATED}' \
                | bcftools filter \
                    --soft-filter LowGQ \
                    -m + \
                    -e 'FMT/GQ < ${INV_MIN_GQ}' \
                | bcftools filter \
                    --soft-filter LowDP \
                    -m + \
                    -e 'FMT/DP < ${INV_MIN_DP}' \
                | bcftools filter \
                    --soft-filter HIGH_DEPTH_GC \
                    -m + \
                    -T '${HIGH_DEPTH_BED}' \
                | bcftools filter \
                    --soft-filter INACCESSIBLE \
                    -m + \
                    -T '${INACCESSIBLE_BED}' \
                    -Oz -o '${SOFTFILTER_OUT}'
            "
        fi
        mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$SOFTFILTER_OUT"
        echo "  [${MODE}] Softfilter VCF written."
    else
        echo "  [${MODE}] Softfilter VCF exists, skipping."
    fi

    echo "  [${MODE}] Producing PASS-only output..."
    if [[ ! -f "$HARDFILT_OUT" ]]; then
        mamba run --no-banner -p "$ENV_BUDDY" \
            bcftools view -f PASS \
                "$SOFTFILTER_OUT" \
                -Oz -o "$HARDFILT_OUT"
        mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$HARDFILT_OUT"
        echo "  [${MODE}] PASS VCF written: ${HARDFILT_OUT}"
    else
        echo "  [${MODE}] PASS VCF exists, skipping."
    fi

    # Companion TSV via SnpSift
    local TSV_OUT="${HARDFILT_OUT%.vcf.gz}.tsv"
    if [[ ! -f "$TSV_OUT" ]]; then
        "$JAVA23" -jar "$SNPSIFT_JAR" extractFields \
            "$HARDFILT_OUT" \
            CHROM POS REF ALT \
            "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].BIOTYPE" \
            INFO/REP \
            > "$TSV_OUT" 2>/dev/null || true
        echo "  [${MODE}] TSV written: ${TSV_OUT}"
    fi
}

filter_and_annotate \
    "$VARIANTS_SNPEFF" \
    "variant" \
    "$OUTPUT_DIR/variants_softfilter.vcf.gz" \
    "$OUTPUT_DIR/variants_filt.vcf.gz"

filter_and_annotate \
    "$INVARIANTS_SNPEFF" \
    "invariant" \
    "$OUTPUT_DIR/invariant_softfilter.vcf.gz" \
    "$OUTPUT_DIR/invariant_filt.vcf.gz"

# ─────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────
echo ""
echo "=== Filter summary for CHR=${CHR} ==="

{
    echo "=== CHR=${CHR} filter summary $(date) ==="
    echo ""
    for VCF_LABEL in variants_filt variants_softfilter invariant_filt invariant_softfilter; do
        VCF="$OUTPUT_DIR/${VCF_LABEL}.vcf.gz"
        echo "--- ${VCF_LABEL} ---"
        mamba run --no-banner -p "$ENV_BUDDY" bcftools stats "$VCF" | grep "^SN"
        mamba run --no-banner -p "$ENV_BUDDY" bcftools query -f '%FILTER\n' "$VCF" \
            | sort | uniq -c | sort -rn
        echo ""
    done
} | tee "$OUTPUT_DIR/qc/filter_summary.txt"

echo ""
echo "All outputs written to ${OUTPUT_DIR}"
