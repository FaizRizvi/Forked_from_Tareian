#!/bin/bash
for fn in /Volumes/NOSTRADOMUS/1000_RNA/ERR*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
kallisto quant -i /Users/caz3so/Dropbox/thesis/data/20180620_Salmon_Kallisto/hg19_index.idx -o ${samp}_quant ${fn}/${samp}_1.fastq.gz ${fn}/${samp}_2.fastq.gz
done