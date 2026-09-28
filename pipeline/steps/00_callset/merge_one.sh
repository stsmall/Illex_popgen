#!/usr/bin/env bash
# GT-ONLY combined all-sites merge for one chromosome.
# Fixes the bcftools-merge corruption at ALT=* / phased-indel sites: the phasing
# FORMAT fields (PGT/PID/PS) get misaligned by merge and crash re-parsing. We
# strip each SOURCE to GT-only (INFO kept: QUAL/MQ/FS survive) BEFORE merge —
# stripping after merge is impossible because the corrupt records crash on read.
# GT-only per user decision (2026-07-02): depth/GL info lives in the ANGSD plane.
# SHARED MACHINE cap: 100 cores (driver runs this at -P12).
set -uo pipefail
c=$1
source /sietch_colab/data_share/illex/popgen_data/analysis/config.sh
BB=/home/ssmall/miniforge3/envs/bioinfo-buddy
cd "$ANA/steps/00_callset"; mkdir -p merged tmp_strip
out=merged/all.$c.vcf.gz
[ -f "$out.tbi" ] && { echo "[$(date)] skip $c (done)"; exit 0; }
rm -f "$out"
bk=tmp_strip/baker.$c.gt.vcf.gz
sm=tmp_strip/small.$c.gt.vcf.gz
echo "[$(date)] START chr $c (GT-only)"
mamba run --no-banner -p "$BB" bash -c "
  set -o pipefail
  bcftools annotate -x '^FORMAT/GT' '$BAKER_VCFDIR/all.$c.vcf.gz' -Oz -o '$bk' && bcftools index -t '$bk'
  bcftools annotate -x '^FORMAT/GT' '$SMALL_VCFDIR/all.$c.vcf.gz' -Oz -o '$sm' && bcftools index -t '$sm'
  bcftools merge -m all '$bk' '$sm' -Ou \
    | bcftools view -S ingroup.txt --force-samples -Oz -o '$out'
  bcftools index -t '$out'
" && { rm -f "$bk" "$bk.tbi" "$sm" "$sm.tbi"; echo "[$(date)] DONE chr $c"; } \
  || { echo "[$(date)] FAILED chr $c"; rm -f "$out" "$bk" "$bk.tbi" "$sm" "$sm.tbi"; exit 1; }
