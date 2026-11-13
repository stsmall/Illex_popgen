#!/usr/bin/env bash
# Inputs for karyotype_qc.py (review M4), run in analysis/steps/03_karyotype/m4_qc.
#  1. chr20 biallelic PASS SNPs of the 619 -> plink2 --hardy/--het/--missing (all SNPs and MAF >= 0.05).
#  2. Depth-matched windows: karyotype_qc.py match, then compute_flank_div_snps.py with the matched groups.
#  3. Genotype likelihoods at the biallelic near-fixed sites between the breakpoints (ANGSD -doGlf 2,
#     -doMajorMinor 3 with REF/ALT from the callset), then karyotype_qc.py glpost.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
M=$(cd "$HERE/../../../.." && pwd)
A=/sietch_colab/data_share/illex/popgen_data/analysis/steps
O=$A/03_karyotype/m4_qc; mkdir -p "$O"; cd "$O"
B=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
P=/home/ssmall/bin/plink2
PY=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/python
ANGSD=/home/ssmall/programs/angsd/angsd
awk '!/^#/ && $1!="sample"{print $1}' "$M/data/qc/karyotypes_619.tsv" > keep_619.txt
# 1
$B view -S keep_619.txt -f PASS -m2 -M2 -v snps --threads 4 -Ob -o chr20_bi.bcf "$A/00_callset/filtered/20/variants_filt.vcf.gz"
$P --bcf chr20_bi.bcf --set-all-var-ids @:# --make-pgen --threads 8 --out chr20 > /dev/null
$P --pfile chr20 --hardy --het --missing --threads 8 --out chr20_all > /dev/null
$P --pfile chr20 --maf 0.05 --hardy --het --missing --threads 8 --out chr20_maf05 > /dev/null
rm -f chr20_bi.bcf
$PY "$HERE/karyotype_qc.py" het | tee het_summary.txt
# 2
$PY "$HERE/karyotype_qc.py" match
POPS=$O/pops_depthmatched.tsv OUT=$O/chr2_40_100_windows_depthmatched.csv CUDA_VISIBLE_DEVICES=${GPU:-1} \
  /home/ssmall/miniforge3/envs/varbuddy-pggpu/bin/python "$M/pipeline/steps/03_karyotype/inversion_content/flank_div/compute_flank_div_snps.py"
$PY "$HERE/karyotype_qc.py" ratios | tee ratios_depthmatched.txt
# 3
$B query -r 2:60540000-79500000 -i 'TYPE="snp" && N_ALT=1' -f '%POS\t%REF\t%ALT\n' \
  "$A/00_callset/filtered/2/variants_filt.vcf.gz" > inv_biallelic.tsv
$PY "$HERE/karyotype_qc.py" sites
grep -f <(sed 's/$/.bam/' keep_619.txt) "$A/04_angsd_chr2/bamlists/bakerMeta_bamlist.txt" > bams619.txt
$ANGSD sites index nearfixed_sites.txt
$ANGSD -b bams619.txt -rf nearfixed_regions.txt -sites nearfixed_sites.txt -doMajorMinor 3 -GL 1 -doGlf 2 \
  -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -nThreads 8 -out gl_nearfixed
$PY "$HERE/karyotype_qc.py" glpost | tee gl_posterior_summary.txt
