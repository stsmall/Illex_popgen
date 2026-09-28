#!/usr/bin/env bash
# filter_variants_one.sh — per-chrom VARIANT filter driver (Task 0.2 re-arch)
#
# Usage:
#   CHR=Z bash filter_variants_one.sh
#   CHR=Z REGION=Z:1-3000000 bash filter_variants_one.sh   # slice mode for testing
#
# Outputs (filtered/<chr>/):
#   variants_softfilter.vcf.gz   — all variants, FILTER column populated
#   variants_filt.vcf.gz         — PASS-only
#   variants_filt.tsv            — SnpSift extractFields TSV
#   qc/exchet_distribution.pdf   — ExcHet QC plot
#
# Why this is separate from the invariant driver:
#   merged/all.<chr>.vcf.gz are all-sites (~99.9% invariant).
#   fill-tags, snpEff, and REP are only meaningful on variants.
#   Running them over the invariant stream causes multi-day runtimes.
#   This driver splits FIRST (cheap bcftools view -c1) then annotates.
#
# Key envs (per interface-notes §7):
#   bioinfo-buddy  — bcftools 1.22 + fill-tags plugin + bedtools + bgzip/tabix
#   JAVA23         — Java 23+ for snpEff 5.4c (mkado-vcf Java 11 is too old)
#
# BUG FIX vs run_filter.sh:
#   Original used -T <bed> in bcftools filter, which is --targets-file and SUBSETS
#   output to BED sites only (drops everything outside the BED).
#   Correct flag for soft-filtering regions: --mask-file <bed>
#   Verified: -T Z:1-10000 test: 101 records in → 54 with -T (subset), 101 with --mask-file.

set -euo pipefail

# ─── Config ────────────────────────────────────────────────────────────────────
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh

[[ -z "${CHR:-}" ]] && { echo "ERROR: CHR must be set (e.g. CHR=Z bash $0)" >&2; exit 1; }

SHARED_DIR="$ANA/steps/00_callset/filtered/_shared"
OUTPUT_DIR="$ANA/steps/00_callset/filtered/${CHR}"
INPUT_VCF="$ANA/steps/00_callset/merged/all.${CHR}.vcf.gz"
REP_BED=/sietch_colab/data_share/illex/annotations/repeat_masking_dir/temp_/repeats.merged.bed

# snpEff
SNPEFF_DB="illex"
SNPEFF_DATA_DIR="/home/ssmall/programs/snpEff/data"
SNPEFF_JAR="/home/ssmall/programs/snpEff/snpEff.jar"
SNPSIFT_JAR="/home/ssmall/programs/snpEff/SnpSift.jar"
JAVA23="/home/ssmall/miniforge3/envs/mkado_env/lib/jvm/bin/java"

ENV_BUDDY=/home/ssmall/miniforge3/envs/bioinfo-buddy
THREADS=16

# Filter thresholds — variants (GT-only callset: no per-sample DP/GQ/AD available,
# and per-sample site-filters flag ~100% of sites in a 653-sample 9x cohort.
# Per 2026-07-02 decision: keep site QUAL/MQ/FS + PARALOG(ExcHet) + region masks
# + site call-rate F_MISSING; DROP LowGQ/LowDP/BadAB.)
VAR_MIN_QUAL=20
VAR_MIN_MQ=40
VAR_MAX_FS=60
VAR_EXCHET_PARALOG=0.001
VAR_MIN_AF_PARALOG=0.1
VAR_MAX_AF_PARALOG=0.9
VAR_MAX_FMISSING=0.5

HIGH_DEPTH_BED="$SHARED_DIR/beds/high_depth_windows.bed"
INACCESSIBLE_BED="$SHARED_DIR/beds/inaccessible.bed"
ANN_HDR="$SHARED_DIR/annotation/extra_info_header.hdr"

# Optional region for slice testing (e.g. REGION=Z:1-3000000)
REGION="${REGION:-}"

# ─── Validate inputs ────────────────────────────────────────────────────────────
for f in "$HIGH_DEPTH_BED" "$INACCESSIBLE_BED" "$ANN_HDR" "$INPUT_VCF" "$REP_BED"; do
    [[ -f "$f" ]] || { echo "ERROR: missing required file: $f" >&2; exit 1; }
done

# ─── Setup dirs ─────────────────────────────────────────────────────────────────
mkdir -p "$OUTPUT_DIR"/{tmp,qc,annotation}

echo "============================================================"
echo "filter_variants_one.sh  CHR=${CHR}${REGION:+  REGION=${REGION}}"
echo "INPUT:  ${INPUT_VCF}"
echo "OUTPUT: ${OUTPUT_DIR}"
echo "============================================================"

# ─── STEP 1: Extract variants and fill tags ─────────────────────────────────────
# Split FIRST so fill-tags only runs on the tiny variant fraction (~0.1%)
VARIANTS_TAGGED="$OUTPUT_DIR/tmp/01_variants_tagged.vcf.gz"

# DATA-INTEGRITY: the merged VCFs contain structurally malformed records (truncated
# FORMAT sample fields) at some ALT=* spanning-deletion sites (e.g. Z:24852210+).
# bcftools view CRASHES parsing these ("vcf_parse_format_fill5: Couldn't read GT data")
# before any -c1/-e/-C0 filter applies. An awk text-sanitizer works but re-serializing
# 653-sample rows is prohibitively slow (~40min/5Mb → hours per chromosome).
# RESOLUTION: fix ONCE upstream. If a sanitized clean file exists it is preferred; the
# drivers then use the FAST indexed bcftools view path. See report §0 for the fix command.
CLEAN_VCF="$ANA/steps/00_callset/merged/all.${CHR}.clean.vcf.gz"
if [[ -f "$CLEAN_VCF" ]]; then
    SRC_VCF="$CLEAN_VCF"
    echo "  Using sanitized input: $CLEAN_VCF"
else
    SRC_VCF="$INPUT_VCF"
fi

echo "[1/6] Extract variants + fill tags (AF, HWE, ExcHet)"
if [[ ! -f "$VARIANTS_TAGGED" ]]; then
    REGION_ARG=""
    [[ -n "$REGION" ]] && REGION_ARG="-r $REGION"
    # Fast indexed path. NOTE: crashes on the raw merged VCF if the corrupt ALT=* records
    # (see above) are within range — provide all.<chr>.clean.vcf.gz to avoid this.
    mamba run --no-banner -p "$ENV_BUDDY" bash -c "
        bcftools view -c1 ${REGION_ARG} '${SRC_VCF}' -Ou \
        | bcftools +fill-tags -Oz -o '${VARIANTS_TAGGED}' -- -t AF,HWE,ExcHet,F_MISSING
    "
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$VARIANTS_TAGGED"
    echo "  Variants tagged: $(mamba run --no-banner -p "$ENV_BUDDY" bcftools stats "$VARIANTS_TAGGED" | awk -F'\t' '/^SN.*records/{print $4}') records"
else
    echo "  Exists, skipping."
fi

# ─── STEP 2: snpEff annotation ──────────────────────────────────────────────────
VARIANTS_SNPEFF="$OUTPUT_DIR/tmp/02_variants_snpeff.vcf.gz"

echo "[2/6] snpEff annotation (DB: ${SNPEFF_DB})"
if [[ ! -f "$VARIANTS_SNPEFF" ]]; then
    # -noCheckCds / -noCheckProtein are BUILD-time flags only; not valid for ann.
    # Write uncompressed VCF first, then use bcftools view -Oz to compress.
    # Cannot pipe java stdout directly to mamba-env bgzip: produces corrupted output
    # (bgzip outside mamba env mangles bytes; bcftools view inside env is reliable).
    SNPEFF_TMP="$OUTPUT_DIR/tmp/02_variants_snpeff.tmp.vcf"
    "$JAVA23" -Xmx16g -jar "$SNPEFF_JAR" ann \
        -dataDir "$SNPEFF_DATA_DIR" \
        -noStats \
        -noLog \
        "$SNPEFF_DB" \
        "$VARIANTS_TAGGED" \
        > "$SNPEFF_TMP"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools view "$SNPEFF_TMP" -Oz -o "$VARIANTS_SNPEFF"
    rm -f "$SNPEFF_TMP"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$VARIANTS_SNPEFF"
    echo "  snpEff done."
else
    echo "  Exists, skipping."
fi

# ─── STEP 3: Build per-chrom REP annotation table ───────────────────────────────
# bedtools intersect VCF positions against repeats.merged.bed
# Output: bgzipped CHROM/POS/REF/ALT/REP table for bcftools annotate
# NOTE: runs on variant positions only (much faster than all-sites)
REPEAT_ANN_TAB="$OUTPUT_DIR/annotation/repeat_annotation.tab.gz"

echo "[3/6] Building REP annotation table for chr ${CHR}"
if [[ ! -f "$REPEAT_ANN_TAB" ]]; then
    mamba run --no-banner -p "$ENV_BUDDY" bash -c "
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' '${VARIANTS_SNPEFF}' \
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
    echo "  REP table built."
else
    echo "  Exists, skipping."
fi

# ─── STEP 4: Annotate REP into VCF ─────────────────────────────────────────────
VARIANTS_REP="$OUTPUT_DIR/tmp/03_variants_rep.vcf.gz"

echo "[4/6] Adding REP annotation"
if [[ ! -f "$VARIANTS_REP" ]]; then
    mamba run --no-banner -p "$ENV_BUDDY" \
        bcftools annotate \
            -a "$REPEAT_ANN_TAB" \
            -c CHROM,POS,REF,ALT,INFO/REP \
            -h "$ANN_HDR" \
            "$VARIANTS_SNPEFF" \
            -Oz -o "$VARIANTS_REP"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$VARIANTS_REP"
    echo "  REP annotation added."
else
    echo "  Exists, skipping."
fi

# ─── STEP 5: Soft-filter chain ──────────────────────────────────────────────────
# CRITICAL: Use --mask-file (not -T/--targets-file) for region-based tagging.
# -T subsets output to BED sites only (WRONG). --mask-file tags and keeps all records.
SOFTFILTER_OUT="$OUTPUT_DIR/variants_softfilter.vcf.gz"

echo "[5/6] Applying soft-filter chain"
if [[ ! -f "$SOFTFILTER_OUT" ]]; then
    mamba run --no-banner -p "$ENV_BUDDY" bash -c "
        bcftools norm -m -any '${VARIANTS_REP}' -Ou \
        | bcftools filter \
            --soft-filter LowSiteQual \
            -m+ \
            -e 'QUAL < ${VAR_MIN_QUAL} || MQ < ${VAR_MIN_MQ} || FS > ${VAR_MAX_FS}' \
            -Ou \
        | bcftools filter \
            --soft-filter PARALOG \
            -m+ \
            -e 'AF > ${VAR_MIN_AF_PARALOG} && AF < ${VAR_MAX_AF_PARALOG} && ExcHet < ${VAR_EXCHET_PARALOG}' \
            -Ou \
        | bcftools filter \
            --soft-filter MISSING \
            -m+ \
            -e 'F_MISSING >= ${VAR_MAX_FMISSING}' \
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
    echo "  Softfilter VCF written."
else
    echo "  Exists, skipping."
fi

# PASS-only hard filter
HARDFILT_OUT="$OUTPUT_DIR/variants_filt.vcf.gz"
if [[ ! -f "$HARDFILT_OUT" ]]; then
    mamba run --no-banner -p "$ENV_BUDDY" \
        bcftools view -f PASS "$SOFTFILTER_OUT" -Oz -o "$HARDFILT_OUT"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools index --tbi "$HARDFILT_OUT"
    echo "  PASS-only VCF: ${HARDFILT_OUT}"
fi

# Companion TSV
TSV_OUT="$OUTPUT_DIR/variants_filt.tsv"
if [[ ! -f "$TSV_OUT" ]]; then
    "$JAVA23" -jar "$SNPSIFT_JAR" extractFields \
        "$HARDFILT_OUT" \
        CHROM POS REF ALT \
        "ANN[0].EFFECT" "ANN[0].GENE" "ANN[0].BIOTYPE" \
        INFO/REP \
        > "$TSV_OUT" 2>/dev/null || true
    echo "  TSV written: ${TSV_OUT}"
fi

# ─── STEP 6: QC — ExcHet distribution ──────────────────────────────────────────
QC_PLOT="$OUTPUT_DIR/qc/exchet_distribution.pdf"

echo "[6/6] ExcHet QC"
if [[ ! -f "$QC_PLOT" ]]; then
    QC_DATA="$OUTPUT_DIR/tmp/exchet_data.tsv"
    mamba run --no-banner -p "$ENV_BUDDY" \
        bcftools query \
            -i "AF > ${VAR_MIN_AF_PARALOG} && AF < ${VAR_MAX_AF_PARALOG}" \
            -f '%CHROM\t%POS\t%AF\t%ExcHet\t%INFO/DP\n' \
            "$VARIANTS_TAGGED" \
        > "$QC_DATA"

    QC_RSCRIPT="$OUTPUT_DIR/tmp/exchet_qc.R"
    cat > "$QC_RSCRIPT" << 'REOF'
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outpdf <- args[2]
d <- read.table(infile, header=FALSE,
                col.names=c("chr","pos","af","exchet","dp"))
# af and exchet may be read as character if bcftools outputs "." for missing
d$af     <- suppressWarnings(as.numeric(d$af))
d$exchet <- suppressWarnings(as.numeric(d$exchet))
d <- d[!is.na(d$af) & !is.na(d$exchet), ]

cat("Mid-AF sites (0.1<AF<0.9):", nrow(d), "\n")
if (nrow(d) == 0) { stop("No mid-AF sites found for QC plot.") }

n_below <- sum(d$exchet < 0.001, na.rm=TRUE)
cat("Sites with ExcHet < 0.001:", n_below, "(", round(100*n_below/nrow(d),2), "%)\n")
cat("ExcHet quantiles (mid-AF sites):\n")
print(quantile(d$exchet, c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), na.rm=TRUE))

pdf(outpdf)
par(mfrow=c(1,2))
hist(-log10(d$exchet + 1e-300), breaks=100,
     main=paste0("ExcHet (mid-AF sites, n=", nrow(d), ")"),
     xlab="-log10(ExcHet)")
abline(v=3, col="red", lty=2)
legend("topright", legend="p=0.001", col="red", lty=2)
plot(d$af, -log10(d$exchet + 1e-300), pch=20, cex=0.3,
     xlab="AF", ylab="-log10(ExcHet)", main="ExcHet vs AF")
abline(h=3, col="red", lty=2)
dev.off()
cat("QC plot written to:", outpdf, "\n")
REOF

    Rscript "$QC_RSCRIPT" "$QC_DATA" "$QC_PLOT"
    echo "  QC plot: $QC_PLOT"
else
    echo "  QC plot exists, skipping."
fi

# ─── Summary ────────────────────────────────────────────────────────────────────
echo ""
echo "=== Variant filter summary CHR=${CHR} $(date) ==="
for label in variants_filt variants_softfilter; do
    VCF="$OUTPUT_DIR/${label}.vcf.gz"
    [[ -f "$VCF" ]] || continue
    echo "--- ${label} ---"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools stats "$VCF" | grep "^SN"
    mamba run --no-banner -p "$ENV_BUDDY" bcftools query -f '%FILTER\n' "$VCF" \
        | sort | uniq -c | sort -rn
    echo ""
done | tee "$OUTPUT_DIR/qc/variant_filter_summary.txt"

echo "Done. Outputs in: ${OUTPUT_DIR}"
