#!/usr/bin/env bash
# Full-genome VARIANT filtering: run filter_variants_one.sh per chrom.
# SHARED MACHINE cap: 100 cores. -P12 x (~1-2 bcftools + 1 snpEff java) ~= ~36 cores.
# Skips chroms whose variants_filt.vcf.gz already exists (e.g. chr24 test).
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
cd "$ANA/steps/00_callset"
run_one() {
  local c=$1
  if [ -f "filtered/$c/variants_filt.vcf.gz" ]; then echo "[$(date)] skip $c (variants_filt exists)"; return 0; fi
  echo "[$(date)] START variant-filter chr $c"
  if CHR=$c bash filter_variants_one.sh > "logs/filter_var_$c.log" 2>&1; then
    echo "[$(date)] DONE chr $c"
  else
    echo "[$(date)] FAILED chr $c (see logs/filter_var_$c.log)"
  fi
}
export -f run_one 2>/dev/null || true
export ANA
echo "[$(date)] variant filtering start (-P12)"
printf '%s\n' $CHR_ALL | xargs -P 12 -I{} bash -c 'run_one "$@"' _ {}
echo "[$(date)] variant filtering driver done"
echo "=== per-chrom PASS counts ==="
for c in $CHR_ALL; do
  f="filtered/$c/variants_filt.vcf.gz"
  [ -f "$f" ] && echo "chr $c: $(mamba run --no-banner -p /home/ssmall/miniforge3/envs/bioinfo-buddy bcftools index -n "$f") PASS" || echo "chr $c: MISSING"
done | tee logs/variant_pass_counts.txt
