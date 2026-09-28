#!/usr/bin/env bash
# Empirical diploSHIC genome scan on the baker source VCFs. Per chrom: fvecVcf (--targetPop illex
# selects the SAME 350 illex the model was trained on, subsets BEFORE stats -> n=700 hap match;
# NO downsample) -> predict (real data, NO --simData). 43 autosomes (ensemble set, excl chr2 inv + chr42).
# Run detached:  setsid bash scan_empirical.sh >/dev/null 2>&1 </dev/null &
set -uo pipefail
B=/sietch_colab/data_share/illex/popgen_data
D=$B/analysis/steps/14_sweep_seqmodel
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
SAM=/home/ssmall/bin/samtools
VCFDIR=$B/seq_data/baker_2025/vcfs
FAI=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa.fai
S2P=$D/results/task17_maskval/s2p.tsv
MASKG=$D/../13_diploshic/mask.fa
MODEL=$D/results/fullrun/illexModel
OUT=$D/results/empirical_scan; MASKD=$OUT/masks
mkdir -p "$OUT"/{fvec,preds,logs} "$MASKD"
CHROMS="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"
CONC=6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 CUDA_VISIBLE_DEVICES=""
export DE SAM VCFDIR FAI S2P MASKG MASKD MODEL OUT
PROG=$OUT/scan.progress
say(){ echo "$(date -u +%FT%TZ) $*" >> "$PROG"; }
say "SCAN_START conc=$CONC"

# per-chrom small mask FASTA (once) so fvecVcf doesn't re-read the 3.3G genome mask each call
for c in $CHROMS; do [ -s "$MASKD/mask.$c.fa" ] || "$SAM" faidx "$MASKG" "$c" > "$MASKD/mask.$c.fa" 2>/dev/null; done
say "MASKS_READY ($(ls "$MASKD"/mask.*.fa 2>/dev/null | wc -l)/43)"

do_chrom(){
  local c=$1
  local vcf=$VCFDIR/all.$c.vcf.gz
  local fvec=$OUT/fvec/chr$c.fvec out=$OUT/preds/chr$c.preds lock=$OUT/preds/chr$c.lock
  [ -s "$out" ] && { echo "skip $c"; return; }
  mkdir "$lock" 2>/dev/null || return
  local L; L=$(awk -v c="$c" '$1==c{print $2}' "$FAI")
  if [ ! -s "$fvec" ]; then
    "$DE/diploSHIC" fvecVcf diploid "$vcf" "$c" "$L" "$fvec" \
      --targetPop illex --sampleToPopFileName "$S2P" \
      --winSize 1100000 --numSubWins 11 --maskFileName "$MASKD/mask.$c.fa" \
      --unmaskedFracCutoff 0.25 --unmaskedGenoFracCutoff 0.5 > "$OUT/logs/chr$c.fvec.log" 2>&1
  fi
  if [ ! -s "$fvec" ]; then echo "FVEC_FAIL $c"; rmdir "$lock" 2>/dev/null; return; fi
  "$DE/diploSHIC" predict "$MODEL.json" "$MODEL.weights.h5" "$fvec" "$out" --numSubWins 11 > "$OUT/logs/chr$c.pred.log" 2>&1
  if [ -s "$out" ]; then echo "DONE $c ($(($(wc -l <"$out")-1)) windows)"; else echo "PRED_FAIL $c"; fi
  rmdir "$lock" 2>/dev/null
}
export -f do_chrom
rm -rf "$OUT"/preds/*.lock 2>/dev/null
say "FVEC_PREDICT_START"
printf '%s\n' $CHROMS | xargs -P "$CONC" -I{} bash -c 'do_chrom "$@"' _ {}
rm -rf "$OUT"/preds/*.lock 2>/dev/null
np=$(ls "$OUT"/preds/chr*.preds 2>/dev/null | wc -l)
say "SCAN_DONE ($np/43 chroms predicted)"
# combine into one genome-wide predictions table
first=$(ls "$OUT"/preds/chr*.preds 2>/dev/null | head -1)
if [ -n "$first" ]; then
  { head -1 "$first"; for f in "$OUT"/preds/chr*.preds; do tail -n +2 "$f"; done; } > "$OUT/genome.preds"
  say "COMBINED $OUT/genome.preds ($(($(wc -l <"$OUT/genome.preds")-1)) windows)"
fi
say "ALL_DONE"
