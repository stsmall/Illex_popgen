#!/bin/bash
# folded MAF spectrum shape: block vs control, chr1 and chr34 (thin every 10th record)
for spec in "1 23000000-31000000 chr1_BLOCK" "1 33000000-41000000 chr1_CONTROL" "34 29900000-44700000 chr34_BLOCK" "34 5000000-20000000 chr34_CONTROL"; do
  set -- $spec
  /home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools query -r "$1:$2" -f '%INFO/AF\n' /sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/$1/variants_filt.vcf.gz 2>/dev/null | awk -v lab="$3" 'NR%10==0{split($1,a,","); af=a[1]+0; m=(af>0.5)?1-af:af; n++; if(m<0.05)r++; else if(m<0.2)mid++; else hi++} END{printf "%-14s SNPs=%7d  MAF<0.05: %.3f  0.05-0.2: %.3f  >=0.2: %.3f  (int/rare ratio %.3f)\n",lab,n,r/n,mid/n,hi/n,hi/r}'
done > /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_sfs.txt 2>&1
echo DONE >> /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/block_sfs.txt
