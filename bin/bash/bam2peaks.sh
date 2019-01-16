#!/bin/bash
cores=4

for i in *hg19_biowardrobe.bam;
do
base_file_name=`basename $i .bam`

filtered_bam=${base_file_name}_filtered.bam
name_sorted=${base_file_name}_namesort.bam
fix_mate=${base_file_name}_fixmate.bam
position_sorted=${base_file_name}_positionsorted.bam
deduped=${base_file_name}_deduped.bam
noMito=${base_file_name}_noMito.bam
sam=${base_file_name}.sam
peakdeck_name=${base_file_name}_peakdeck.bed
peakdeck_tmp=${base_file_name}_peaktemp.bed
chrom_sizes=hg19.chrom.sizes.txt

###BODY####

echo "Filtering: In Progress"

# View the BAM file and filter based on the quality. 
samtools view -@ ${cores} -b -q 30 -o ${filtered_bam} ${i}

echo "Filtering: Done"
echo "Sorting: In Progress"

#  Sort the file based on the name
samtools sort -n -o ${name_sorted} ${filtered_bam}

echo "Sorting: Done"
echo "Fixing Mates: In Progress"

# Fix mates
samtools fixmate -m ${name_sorted} ${fix_mate}

echo "Fixing Mates: Done"
echo "Sorting Fixed Mates: In Progress"

# Sort files by position now
samtools sort -o ${position_sorted} ${fix_mate}

echo "Sorting Fixed Mates: Done"
echo "Marking Duplicated: In Progress"

# Mark duplicates
samtools markdup -r -s ${position_sorted} ${deduped}

echo "Marking Duplicated: Done"
echo "Indexing Marked Duplicates File: In Progress"

# Index duplicate file
samtools index ${deduped}

echo "Indexing Marked Duplicates File: Done"
echo "Remove Mitochondrial Reads: In Progress"

#samtools idxstats ${deduped} | cut -f 1 | grep -v MT | xargs samtools view -b ${deduped} > ${noMito}
#bam filter ${deduped}.bam ${noMito}.bam -excluderef chrM 
samtools view -@ ${cores} -b ${deduped} chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr1 chr20 chr21 chr22 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY > ${noMito}

echo "Remove Mitochondrial Reads: Done"
echo "Convert BAM to SAM: In Progress"

# Convert to SAM
samtools view -@ ${cores} -h -o ${sam} ${noMito}

echo "Convert BAM to SAM: Done"
echo "Run PeakDeck: In Progress"

# Run PeakDeck
perl /Users/caz3so/workspaces/thesis/bin/perl/peakdeck.pl -P ${sam} -g ${chrom_sizes} -PVAL ON -bin 75 -STEP 25 -back 10000 -npBack 100000 > ${peakdeck_tmp}
perl /Users/caz3so/workspaces/thesis/bin/perl/peakdeck.pl -P SRR891269_hg19_biowardrobe.sam -g hg19.chrom.sizes.txt -PVAL ON -bin 75 -STEP 25 -back 10000 -npBack 100000 > SRR891269_hg19_biowardrobe_temp.bed

awk 'NF == 6{print $1"\t"gensub(/\.5$/, "", 1, $2)"\t"gensub(/\.5$/, "", 1, $3)"\t"$4"\t"$5"\t"$6}' ${peakdeck_tmp} > ${peakdeck_name}

echo "Run PeakDeck: Done"
echo "Run MACS2: In Progress"

# Run MACS2
macs2 callpeak -t ${noMito} -f BAMPE -g hs --nomodel --shift -100 --extsize 200

echo "Run MACS2: Done"
echo "BAM2PEAKS DONE!"
done