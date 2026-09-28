#!/usr/bin/env bash
# Empirical diploSHIC scan with the fvecVcf-RETRAINED model (illexModel_vcf), matched to the
# training ascertainment: per chrom subset illex350 + MAF>0.01 + strip to GT-only (fast) ->
# fvecVcf (real per-chrom mask, same flags as training ms_to_vcf_fvec) -> predict. Produces
# per-window class probs; downstream uses the sweep score (p_hard+p_soft) as an empirical-outlier
# scan (option 1). Reuses the per-chrom mask FASTAs from the old scan.
# Run detached:  setsid bash run_empirical_scan_maf.sh >/dev/null 2>&1 </dev/null &
set -uo pipefail
B=/sietch_colab/data_share/illex/popgen_data
D=$B/analysis/steps/14_sweep_seqmodel
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
BCF=/home/ssmall/bin/bcftools; SAM=/home/ssmall/bin/samtools
VCFDIR=$B/seq_data/baker_2025/vcfs
FAI=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa.fai
S2P=$D/results/task17_maskval/s2p.tsv
MODEL=$D/results/vcfretrain_full/illexModel_vcf
MASKD=$D/results/empirical_scan/masks             # per-chrom mask FASTAs (already built)
OUT=$D/results/empirical_scan_maf; mkdir -p "$OUT"/{vcf,fvec,preds,logs}
CHROMS="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"
CONC=8
export DE BCF VCFDIR FAI S2P MODEL MASKD OUT
ILLEX=$OUT/illex350.txt; awk '$2=="illex"{print $1}' "$S2P" > "$ILLEX"; export ILLEX
PROG=$OUT/scan.progress; : > "$PROG"; say(){ echo "$(date -u +%FT%TZ) $*" >> "$PROG"; }
say "SCANMAF_START conc=$CONC model=$MODEL"

do_chrom(){
  local c=$1
  local vcf=$VCFDIR/all.$c.vcf.gz maf=$OUT/vcf/chr$c.maf.vcf.gz
  local fvec=$OUT/fvec/chr$c.fvec out=$OUT/preds/chr$c.preds lock=$OUT/preds/chr$c.lock
  [ -s "$out" ] && { echo "skip $c"; return; }
  mkdir "$lock" 2>/dev/null || return
  local L; L=$(awk -v c="$c" '$1==c{print $2}' "$FAI")
  if [ ! -s "$maf" ]; then
    "$BCF" view -S "$ILLEX" "$vcf" -Ou 2>/dev/null \
      | "$BCF" view -q 0.01:minor -m2 -M2 -v snps -Ou 2>/dev/null \
      | "$BCF" annotate -x '^FORMAT/GT,INFO' -Oz -o "$maf" 2>/dev/null
    "$BCF" index -t "$maf" 2>/dev/null
  fi
  if [ ! -s "$fvec" ]; then
    "$DE/diploSHIC" fvecVcf diploid "$maf" "$c" "$L" "$fvec" \
      --targetPop illex --sampleToPopFileName "$S2P" \
      --winSize 1100000 --numSubWins 11 --maskFileName "$MASKD/mask.$c.fa" \
      --unmaskedFracCutoff 0.25 --unmaskedGenoFracCutoff 0.5 > "$OUT/logs/chr$c.fvec.log" 2>&1
  fi
  if [ ! -s "$fvec" ]; then echo "FVEC_FAIL $c"; rmdir "$lock" 2>/dev/null; return; fi
  "$DE/diploSHIC" predict "$MODEL.json" "$MODEL.weights.h5" "$fvec" "$out" --numSubWins 11 > "$OUT/logs/chr$c.pred.log" 2>&1
  [ -s "$out" ] && echo "DONE $c ($(($(wc -l <"$out")-1)) win)" || echo "PRED_FAIL $c"
  rmdir "$lock" 2>/dev/null
}
export -f do_chrom
rm -rf "$OUT"/preds/*.lock 2>/dev/null
say "SCAN_RUN"
printf '%s\n' $CHROMS | xargs -P "$CONC" -I{} bash -c 'do_chrom "$@"' _ {} >> "$PROG" 2>&1
np=$(ls "$OUT"/preds/chr*.preds 2>/dev/null | wc -l)
say "SCAN_DONE ($np/43)"
first=$(ls "$OUT"/preds/chr*.preds 2>/dev/null | head -1)
if [ -n "$first" ]; then
  { head -1 "$first"; for f in "$OUT"/preds/chr*.preds; do tail -n +2 "$f"; done; } > "$OUT/genome.preds"
  say "COMBINED genome.preds ($(($(wc -l <"$OUT/genome.preds")-1)) windows)"
fi
say "ALL_DONE"
