#!/usr/bin/env bash
# Per-chrom on the baker SOURCE VCF: normalize + true-biallelic SNPs (keep FULL SFS, no MAF
# filter) -> unfolded SF2 FreqFile (+ deterministic n0 projection) + polarized VCF for RAiSD.
set -uo pipefail
c=$1
source "$(dirname "$0")/env.sh"
RAW=$VCFDIR/all.$c.vcf.gz
[ -f "$RAW" ] || { echo "POLAR_FAIL $c (no vcf $RAW)"; exit 1; }

# Stage A: normalize (traps padded SNPs/disguised multiallelics), true-biallelic SNPs,
# recompute AC/AN/F_MISSING, keep >=50% call rate. NO MAF filter -> SFS tails retained.
"$BCF" norm -f "$REF" -m +snps "$RAW" -Ou 2>/dev/null \
 | "$BCF" view -m2 -M2 -v snps -Ou 2>/dev/null \
 | "$BCF" +fill-tags -Ou -- -t AC,AN,F_MISSING 2>/dev/null \
 | "$BCF" view -e 'F_MISSING>=0.5' -Oz -o "$R/polarized/$c.snps.vcf.gz" 2>/dev/null
"$BCF" index -t -f "$R/polarized/$c.snps.vcf.gz" 2>/dev/null

# Stage B: annotate ancestral -> unfolded FreqFile (position, x=derived, n, folded=0)
"$BCF" +fill-from-fasta "$R/polarized/$c.snps.vcf.gz" -- -c AA -f "$ANC" -h "$S/aa_header.txt" 2>/dev/null \
 | "$BCF" query -f '%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AA\n' \
 | awk -F'\t' -f "$S/polarize.awk" > "$R/freq/$c.sf2.freq" 2> "$R/freq/$c.polar.stats"

# polarized VCF (REF=ancestral => ALT=derived) for RAiSD's unfolded muSFS
"$BCF" +fill-from-fasta "$R/polarized/$c.snps.vcf.gz" -- -c AA -f "$ANC" -h "$S/aa_header.txt" 2>/dev/null \
 | "$BCF" view -e 'INFO/AA="N" || (INFO/AA!=REF && INFO/AA!=ALT)' -Ou 2>/dev/null \
 | "$BCF" norm -f "$ANC" -c s - -Oz -o "$R/polarized/$c.polarized.vcf.gz" 2>/dev/null
"$BCF" index -t -f "$R/polarized/$c.polarized.vcf.gz" 2>/dev/null

# deterministic projection to constant n0 for SweepFinder2 (random draw injects noise)
"$PY" "$S/project_freq_det.py" "$R/freq/$c.sf2.freq" "$R/freq/$c.n$N0.freq" "$N0" 2>/dev/null

nsnp=$(($(wc -l < "$R/freq/$c.sf2.freq") - 1))
echo "POLAR_OK $c snps=$nsnp $(cat "$R/freq/$c.polar.stats" 2>/dev/null)"
