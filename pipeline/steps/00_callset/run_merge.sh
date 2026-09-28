#!/usr/bin/env bash
# Combined illecebrosus all-sites callset: per-chrom merge of baker + small_2026,
# subset to the 653-sample ingroup. The full merge|view|index pipe runs inside a
# SINGLE `mamba run` (piping a BCF stream between two separate `mamba run`
# invocations fails: "Failed to read from standard input: unknown file type").
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
cd "$ANA/steps/00_callset"; mkdir -p merged
for c in $CHR_ALL; do
  out=merged/all.$c.vcf.gz
  [ -f "$out.tbi" ] && { echo "skip $c"; continue; }
  echo "[$(date)] merge chr $c"
  mamba run --no-banner -p "$ENV_BCF" bash -c "
    set -o pipefail
    bcftools merge -m all -Ou '$BAKER_VCFDIR/all.$c.vcf.gz' '$SMALL_VCFDIR/all.$c.vcf.gz' \
      | bcftools view -S ingroup.txt --force-samples -Oz -o '$out' \
      && bcftools index --tbi '$out'
  " || { echo "[$(date)] FAILED chr $c"; exit 1; }
done
echo "[$(date)] merge done"
