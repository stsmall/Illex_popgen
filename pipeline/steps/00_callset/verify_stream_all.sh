#!/usr/bin/env bash
# Genome-wide parse-check: fully stream every merged chrom to /dev/null and record
# any nonzero exit / parse-error stderr. Confirms no ALT=* crash remains anywhere.
# Capped parallelism (SHARED MACHINE <=100 cores).
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy
D=$ANA/steps/00_callset
cd "$D"; mkdir -p logs/verify
check_one() {
  local c=$1
  mamba run --no-banner -p "$BB" bash -c "bcftools view -H '$D/merged/all.$c.vcf.gz' >/dev/null 2>'$D/logs/verify/$c.err'"
  local rc=$?
  local errs=$(wc -l < "$D/logs/verify/$c.err")
  echo "chr $c: exit=$rc err_lines=$errs"
}
export -f check_one 2>/dev/null || true
export ANA BB D
printf '%s\n' $CHR_ALL | xargs -P 12 -I{} bash -c 'check_one "$@"' _ {} | tee logs/verify/summary.txt
echo "=== VERDICT ==="
if grep -qvE "exit=0 err_lines=0" logs/verify/summary.txt; then
  echo "PROBLEM: some chrom failed clean parse:"; grep -vE "exit=0 err_lines=0" logs/verify/summary.txt
else
  echo "ALL 46 CHROMS PARSE CLEAN (exit=0, 0 errors)"
fi
