#!/usr/bin/env bash
# Called by fvec_real.sh via xargs: fvec_one.sh <k> <chrom>
# Subsets VCF, re-MAF-filters, runs fvecVcf for one subsample x one chromosome.
set -uo pipefail
source /sietch_colab/data_share/illex/popgen_data/analysis/steps/13_diploshic/env.sh

k=$1
c=$2
kk=$(printf '%02d' $k)
d=$WD/fvec/real/sub${kk}
mkdir -p $d

# Get chromosome length from reference fai
len=$(awk -v c="$c" '$1==c{print $2}' ${REF}.fai)
if [ -z "$len" ]; then
  echo "ERROR: chromosome $c not found in ${REF}.fai" >&2
  exit 1
fi

sub=$WD/combined/sub${kk}.$c.vcf.gz

# Subset to the 350, recompute + re-filter MAF>0.01 to match sim ascertainment
$BB/bcftools view -S $WD/subsamples/sub_${kk}.txt $WD/combined/all.${c}.vcf.gz -Ou 2>/dev/null \
  | $BB/bcftools +fill-tags -Ou -- -t MAF 2>/dev/null \
  | $BB/bcftools view -e 'INFO/MAF<0.01' -Oz -o $sub 2>/dev/null
$BB/tabix -f -p vcf $sub

# Build sample-to-pop file for this subsample x chrom
printf '%s\tillex\n' $(cat $WD/subsamples/sub_${kk}.txt) > $d/s2p.${c}.tsv

# Run fvecVcf (note: fvecVcf does NOT accept --totalPhysLen, unlike fvecSim)
$DSHIC fvecVcf diploid $sub $c $len $d/${c}.fvec \
  --numSubWins 11 --winSize 44000 \
  --sampleToPopFileName $d/s2p.${c}.tsv --targetPop illex \
  --maskFileName $WD/mask.fa --unmaskedFracCutoff 0.25
# surface silent failures (the chr3-9 mask bug produced 0-output crashes hidden by 2>/dev/null)
[ -s "$d/${c}.fvec" ] && [ "$(wc -l < "$d/${c}.fvec")" -gt 1 ] || echo "WARN: empty/failed fvec sub${kk} chr${c}" >&2

# Clean up temporary VCF
rm -f $sub ${sub}.tbi $d/s2p.${c}.tsv
