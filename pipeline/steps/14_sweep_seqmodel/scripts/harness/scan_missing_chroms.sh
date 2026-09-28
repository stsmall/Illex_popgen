#!/usr/bin/env bash
# Scan the chromosomes the stale hardcoded CHROMS list omitted from the full-SFS scan.
# chr2 (inversion) + chr42 (AUTOSOMAL -- wrongly excluded on the corrected-away "sex chr42" belief).
# EXACT copy of run_empirical_scan_fullsfs.sh::do_chrom (same fvecVcf flags / mask / cutoffs), just a
# different CHROMS set, writing into the SAME empirical_scan_fullsfs/ dir so genome.preds can be rebuilt
# to include them. chrZ is NOT here -- it is the sex chr (ZW females Z-haploid) and needs a ploidy-aware
# (males-only ZZ) run, decided separately.
# Run detached:  setsid bash scan_missing_chroms.sh >/dev/null 2>&1 </dev/null &
set -uo pipefail
B=/sietch_colab/data_share/illex/popgen_data
D=$B/analysis/steps/14_sweep_seqmodel
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
BCF=/home/ssmall/bin/bcftools; SAM=/home/ssmall/bin/samtools
VCFDIR=$B/seq_data/baker_2025/vcfs
FAI=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa.fai
S2P=$D/results/task17_maskval/s2p.tsv
MASKD=$D/results/empirical_scan/masks
OUT=$D/results/empirical_scan_fullsfs; mkdir -p "$OUT"/{vcf,fvec,logs}
CHROMS="2 42"
CONC=2
export DE BCF VCFDIR FAI S2P MASKD OUT
ILLEX=$OUT/illex350.txt; [ -s "$ILLEX" ] || awk '$2=="illex"{print $1}' "$S2P" > "$ILLEX"; export ILLEX
PROG=$OUT/scan_missing.progress; : > "$PROG"; say(){ echo "$(date -u +%FT%TZ) $*" >> "$PROG"; }
say "MISSING_FVEC_START chroms=[$CHROMS] conc=$CONC"

do_chrom(){
  local c=$1
  local vcf=$VCFDIR/all.$c.vcf.gz nomaf=$OUT/vcf/chr$c.nomaf.vcf.gz
  local fvec=$OUT/fvec/chr$c.fvec lock=$OUT/fvec/chr$c.lock
  [ -s "$fvec" ] && { echo "skip $c"; return; }
  mkdir "$lock" 2>/dev/null || return
  local L; L=$(awk -v c="$c" '$1==c{print $2}' "$FAI")
  if [ ! -s "$nomaf" ]; then
    "$BCF" view -S "$ILLEX" "$vcf" -Ou 2>/dev/null \
      | "$BCF" view -m2 -M2 -v snps -Ou 2>/dev/null \
      | "$BCF" annotate -x '^FORMAT/GT,INFO' -Oz -o "$nomaf" 2>/dev/null
    "$BCF" index -t "$nomaf" 2>/dev/null
  fi
  if [ ! -s "$fvec" ]; then
    "$DE/diploSHIC" fvecVcf diploid "$nomaf" "$c" "$L" "$fvec" \
      --targetPop illex --sampleToPopFileName "$S2P" \
      --winSize 1100000 --numSubWins 11 --maskFileName "$MASKD/mask.$c.fa" \
      --unmaskedFracCutoff 0.25 --unmaskedGenoFracCutoff 0.5 > "$OUT/logs/chr$c.fvec.log" 2>&1
  fi
  if [ -s "$fvec" ]; then echo "DONE $c ($(($(wc -l <"$fvec")-1)) win)"; else echo "FVEC_FAIL $c"; fi
  rmdir "$lock" 2>/dev/null
}
export -f do_chrom
rm -rf "$OUT"/fvec/chr2.lock "$OUT"/fvec/chr42.lock 2>/dev/null
say "FVEC_RUN"
printf '%s\n' $CHROMS | xargs -P "$CONC" -I{} bash -c 'do_chrom "$@"' _ {} >> "$PROG" 2>&1
say "MISSING_FVEC_DONE ($(ls "$OUT"/fvec/chr2.fvec "$OUT"/fvec/chr42.fvec 2>/dev/null | wc -l)/2)"
