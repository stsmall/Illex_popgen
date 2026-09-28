#!/usr/bin/env bash
# Empirical diploSHIC fvec, FULL-SFS (NO MAF filter) -- for the full-SFS retrain (2026-08-11).
# Identical to run_empirical_scan_maf.sh EXCEPT:
#   * drops `bcftools view -q 0.01:minor` (keeps the -m2 -M2 -v snps biallelic restriction)
#   * writes to results/empirical_scan_fullsfs/ (old MAF scan left intact for comparison)
#   * GENERATES FVEC ONLY -- predict is DEFERRED until the full-SFS model is trained.
#     (Predicting full-SFS features with the MAF-trained illexModel_vcf would just recreate
#      the train/scan mismatch, so predict waits for the retrained model.)
# Reuses the per-chrom accessibility mask FASTAs from the old scan. fvecVcf flags identical to
# the training path (ms_to_vcf_fvec.py): winSize 1.1Mb, 11 subwins, unmasked cutoffs 0.25/0.5.
# Run detached:  setsid bash run_empirical_scan_fullsfs.sh >/dev/null 2>&1 </dev/null &
set -uo pipefail
B=/sietch_colab/data_share/illex/popgen_data
D=$B/analysis/steps/14_sweep_seqmodel
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
BCF=/home/ssmall/bin/bcftools; SAM=/home/ssmall/bin/samtools
VCFDIR=$B/seq_data/baker_2025/vcfs
FAI=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa.fai
S2P=$D/results/task17_maskval/s2p_clean350.tsv
MASKD=$D/results/empirical_scan/masks             # per-chrom mask FASTAs (already built)
OUT=$D/results/empirical_scan_clean350; mkdir -p "$OUT"/{vcf,fvec,logs}
CHROMS="1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 43 44 45"
CONC=12
export DE BCF VCFDIR FAI S2P MASKD OUT
ILLEX=$OUT/illex350.txt; awk '$2=="illex"{print $1}' "$S2P" > "$ILLEX"; export ILLEX
PROG=$OUT/scan.progress; : > "$PROG"; say(){ echo "$(date -u +%FT%TZ) $*" >> "$PROG"; }
say "FULLSFS_FVEC_START conc=$CONC (NO MAF filter; fvec only, predict deferred)"

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
rm -rf "$OUT"/fvec/*.lock 2>/dev/null
say "FVEC_RUN"
printf '%s\n' $CHROMS | xargs -P "$CONC" -I{} bash -c 'do_chrom "$@"' _ {} >> "$PROG" 2>&1
nf=$(ls "$OUT"/fvec/chr*.fvec 2>/dev/null | wc -l)
say "FULLSFS_FVEC_DONE ($nf/43 fvecs; predict deferred to full-SFS model)"
