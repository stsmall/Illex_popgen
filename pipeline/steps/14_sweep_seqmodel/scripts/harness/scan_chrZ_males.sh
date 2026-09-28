#!/usr/bin/env bash
# chrZ sweep scan, MALES-ONLY (ZZ = Z-diploid). ZW females are Z-haploid and are dropped entirely
# (a diploid diploSHIC scan cannot use haploid genotypes). Subset = illex ∩ males = 171 diploid samples
# (n=171 vs the model's training n=350 -> SFS-feature baselines shift; state this caveat in methods).
# Mirrors run_empirical_scan_fullsfs.sh::do_chrom EXACTLY except -S uses the Zmales list, and writes
# chrZ.* into the canonical empirical_scan_fullsfs/ dir. Predict is a separate step (full-SFS model).
# Run detached:  setsid bash scan_chrZ_males.sh >/dev/null 2>&1 </dev/null &
set -uo pipefail
B=/sietch_colab/data_share/illex/popgen_data
D=$B/analysis/steps/14_sweep_seqmodel
DE=/home/ssmall/miniforge3/envs/diploshic_env/bin
BCF=/home/ssmall/bin/bcftools
VCFDIR=$B/seq_data/baker_2025/vcfs
FAI=/sietch_colab/data_share/illex/assembly/current_assembly_version/Illex_F24.primary.clean.fa.fai
S2P=$D/results/task17_maskval/s2p.tsv
MASKD=$D/results/empirical_scan/masks
OUT=$D/results/empirical_scan_fullsfs; mkdir -p "$OUT"/{vcf,fvec,logs}
SUB=$OUT/Zmales.illex.list       # 171 illex males (built by hand)
c=Z
PROG=$OUT/scan_chrZ.progress; : > "$PROG"; say(){ echo "$(date -u +%FT%TZ) $*" >> "$PROG"; }
say "CHRZ_MALES_FVEC_START n=$(wc -l < "$SUB")"

vcf=$VCFDIR/all.$c.vcf.gz nomaf=$OUT/vcf/chr$c.nomaf.vcf.gz fvec=$OUT/fvec/chr$c.fvec
L=$(awk -v c="$c" '$1==c{print $2}' "$FAI")
if [ ! -s "$nomaf" ]; then
  say "SUBSET_VCF"
  "$BCF" view -S "$SUB" "$vcf" -Ou 2>>"$OUT/logs/chr$c.subset.log" \
    | "$BCF" view -m2 -M2 -v snps -Ou 2>>"$OUT/logs/chr$c.subset.log" \
    | "$BCF" annotate -x '^FORMAT/GT,INFO' -Oz -o "$nomaf" 2>>"$OUT/logs/chr$c.subset.log"
  "$BCF" index -t "$nomaf" 2>/dev/null
fi
if [ ! -s "$fvec" ]; then
  say "FVEC"
  "$DE/diploSHIC" fvecVcf diploid "$nomaf" "$c" "$L" "$fvec" \
    --targetPop illex --sampleToPopFileName "$S2P" \
    --winSize 1100000 --numSubWins 11 --maskFileName "$MASKD/mask.$c.fa" \
    --unmaskedFracCutoff 0.25 --unmaskedGenoFracCutoff 0.5 > "$OUT/logs/chr$c.fvec.log" 2>&1
fi
if [ -s "$fvec" ]; then say "CHRZ_FVEC_DONE ($(($(wc -l <"$fvec")-1)) win)"; else say "CHRZ_FVEC_FAIL"; fi
