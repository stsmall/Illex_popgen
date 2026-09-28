#!/bin/bash
# het fraction + site depth for chr1 block (23-31 Mb) vs control (33-41 Mb); thin to every 25th record
for reg in "1:23000000-31000000 BLOCK" "1:33000000-41000000 CONTROL"; do
  set -- $reg
  /home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools query -r "$1" -f '%POS\t%INFO/DP[\t%GT]\n' /sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered/1/variants_filt.vcf.gz 2>/dev/null | awk -v lab="$2" 'NR%25==0{
    het=0;hom=0;miss=0; for(i=3;i<=NF;i++){g=$i; if(g=="./."||g==".|."||g==".")miss++; else if(g=="0/1"||g=="1/0"||g=="0|1"||g=="1|0")het++; else hom++}
    called=het+hom; if(called>0){hf=het/called; sum+=hf; n++; if(hf>0.6)hi++; dp+=$2}
  } END{printf "%s: sites=%d  mean het-frac=%.3f  frac sites with >60%% het=%.3f  mean INFO/DP=%.0f\n",lab,n,sum/n,hi/n,dp/n}'
done > /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/chr1_block_het.txt 2>&1
echo DONE >> /sietch_colab/data_share/illex/popgen_data/analysis/steps/03_karyotype/inversion_content/chr1_block_het.txt
