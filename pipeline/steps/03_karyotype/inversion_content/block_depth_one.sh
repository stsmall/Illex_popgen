#!/bin/bash
b=$1; s=$(basename $b .bam)
c1=$(/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/samtools view -c -q 20 -F 0x904 $b 1:23000000-31000000 2>/dev/null); k1=$(/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/samtools view -c -q 20 -F 0x904 $b 1:33000000-41000000 2>/dev/null)
c34=$(/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/samtools view -c -q 20 -F 0x904 $b 34:29900000-44700000 2>/dev/null); k34=$(/home/ssmall/miniforge3/envs/bioinfo-buddy/bin/samtools view -c -q 20 -F 0x904 $b 34:5000000-20000000 2>/dev/null)
echo -e "$s\t$c1\t$k1\t$c34\t$k34"
