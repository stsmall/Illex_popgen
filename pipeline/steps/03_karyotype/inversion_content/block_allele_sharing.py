import subprocess, numpy as np
BCF="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/bcftools"; SAM="/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/samtools"
A="/sietch_colab/data_share/illex/popgen_data/analysis/steps/00_callset/filtered"
ARG="/sietch_colab/data_share/illex/popgen_data/polarize/argentinus_input/argentinus_ref.fa"
COIN="/sietch_colab/data_share/illex/popgen_data/mkado_illex/inputs/coindetii.vcf.gz"
def run(cmd): return subprocess.run(cmd,capture_output=True,text=True).stdout
def region(c,lo,hi,thin=20):
    snps=[l.split("\t") for i,l in enumerate(run([BCF,"query","-r",f"{c}:{lo}-{hi}","-f","%POS\t%REF\t%ALT\t%INFO/AF\n",f"{A}/{c}/variants_filt.vcf.gz"]).splitlines()) if i%thin==0]
    snps=[(int(p),r,a.split(",")[0],float(af.split(",")[0])) for p,r,a,af in snps if len(r)==1 and len(a.split(",")[0])==1]
    fa="".join(run([SAM,"faidx",ARG,f"{c}:{lo}-{hi}"]).split("\n")[1:]).upper()
    coin={int(p):a for p,a in (l.split("\t") for l in run([BCF,"query","-r",f"{c}:{lo}-{hi}","-f","%POS\t%ALT\n",COIN]).splitlines())}
    n=len(snps); argN=argA=argR=0; coinA=coinR=coinO=0; argA_common=0; n_common=0
    for p,r,a,af in snps:
        b=fa[p-lo] if 0<=p-lo<len(fa) else "N"
        if b in "N-": argN+=1
        elif b==a: argA+=1
        elif b==r: argR+=1
        cb=coin.get(p)
        if cb is None: coinR+=1
        elif cb==a: coinA+=1
        else: coinO+=1
        maf=min(af,1-af)
        if maf>=0.2:
            n_common+=1
            if b==a: argA_common+=1
    called=argA+argR
    return dict(n=n,arg_alt=argA/max(called,1),arg_called=called/n,coin_alt=coinA/n,coin_other=coinO/n,arg_alt_common=argA_common/max(n_common,1),n_common=n_common)
for lab,c,lo,hi in (("chr1 BLOCK","1",23000000,31000000),("chr1 CONTROL","1",33000000,41000000),
                    ("chr34 BLOCK","34",29900000,44700000),("chr34 CONTROL","34",5000000,20000000)):
    r=region(c,lo,hi)
    print(f"{lab:13s} SNPs={r['n']:6d} | ALT==argentinus: {100*r['arg_alt']:.1f}% (arg base called at {100*r['arg_called']:.0f}% of SNPs) | among MAF>=0.2 SNPs (n={r['n_common']}): ALT==arg {100*r['arg_alt_common']:.1f}% | ALT==coindetii: {100*r['coin_alt']:.1f}% (coin=other allele {100*r['coin_other']:.1f}%)",flush=True)
print("DONE_SHARING")
