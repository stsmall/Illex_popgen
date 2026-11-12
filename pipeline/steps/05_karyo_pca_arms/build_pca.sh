#!/usr/bin/env bash
# PCA for main-text Fig 2a (autosomes, chromosome 2 excluded) and Fig 2b (chrZ) on the 619 analysed individuals.
# Input: the per-chromosome LD-thinned VCFs of thin_for_pca.sh (baker-633, biallelic SNPs, MAF >= 0.05, 1 SNP per
# 10 kb). The excluded individuals (data/qc/excluded_14.tsv) are removed and MAF >= 0.05 is re-applied on the 619.
# Proportion of variance explained = eigenvalue / trace of the full relationship matrix (written by plink2
# --make-rel), NOT eigenvalue / sum of the top 10 eigenvalues.
# Writes data/{auto,z}_pca.eigenvec, .eigenval and _summary.tsv (n_snps, pve_pc1, pve_pc2).
set -euo pipefail
ANA=/sietch_colab/data_share/illex/popgen_data/analysis
TH=$ANA/steps/05_karyo_pca_arms/thinned
M=$(cd "$(dirname "$0")/../../.." && pwd)
B=/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools
PY=${PY:-python3}
W=$(mktemp -d); trap 'command rm -rf "$W"' EXIT
grep -v '^#' "$M/data/qc/excluded_14.tsv" | cut -f1 | grep -v "^sample$" > "$W/empty.txt"
grep -vxFf "$W/empty.txt" "$ANA/steps/00_callset/baker_633.txt" | awk '{print $1"\t"$1}' > "$W/keep.fid"
$B concat -Oz -o "$W/auto.vcf.gz" $(for c in 1 $(seq 3 45); do echo "$TH/chr$c.thin.vcf.gz"; done) 2>/dev/null
for x in "auto:$W/auto.vcf.gz" "z:$TH/chrZ.thin.vcf.gz"; do
  n=${x%%:*}; v=${x#*:}
  plink2 --vcf "$v" --double-id --chr-set 45 --allow-extra-chr --keep "$W/keep.fid" --maf 0.05 \
    --set-all-var-ids '@:#' --threads 8 --memory 8000 --pca 10 --make-rel square --out "$W/$n" > /dev/null
  cp "$W/$n.eigenvec" "$M/data/${n}_pca.eigenvec"; cp "$W/$n.eigenval" "$M/data/${n}_pca.eigenval"
  nsnp=$(grep -o '[0-9]* variants remaining after main filters' "$W/$n.log" | cut -d' ' -f1)
  $PY - "$W/$n.rel" "$W/$n.eigenval" "$nsnp" > "$M/data/${n}_pca_summary.tsv" <<'PYEOF'
import sys, numpy as np
R = np.loadtxt(sys.argv[1]); ev = np.loadtxt(sys.argv[2]); tr = np.trace(R)
print("n_snps\tpve_pc1\tpve_pc2\ttrace"); print(f"{sys.argv[3]}\t{100*ev[0]/tr:.3f}\t{100*ev[1]/tr:.3f}\t{tr:.2f}")
PYEOF
  cat "$M/data/${n}_pca_summary.tsv"
done
