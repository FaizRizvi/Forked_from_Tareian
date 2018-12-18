#!/bin/bash
for fn in /Volumes/NOSTRADOMUS/1000_RNA/*;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i hg19_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 -o quants/${samp}_quant
done
