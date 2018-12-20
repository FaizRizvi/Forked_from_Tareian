#!/bin/bash

#########Input#########
sample=$1
memory=$2

inDirRoot=/data/miraldiLab/reference_data/Buenrostro/fastq
outDirRoot=/data/miraldiLab/Tareian/results
tmpDirRoot=/data/miraldiLab/Temp
bowtieindex="/data/miraldiLab/reference_data/genome/hg19_index/hg19"

#########Set Parameters from input#########
outDir=${outDirRoot}/${sample}
mkdir -p ${outDir} || { echo "Cannot create \"${outDir}\". Exiting." ; exit 1 ; }
wwi=$(pwd)
wd=$(mktemp -d ${tmpDirRoot}/${USER}_qc_bfsw_XXXXXXXXXX)

fq1=${inDirRoot}/${sample}_1.fastq
fq2=${inDirRoot}/${sample}_2.fastq

sam=${sample}.sam
bam=${sample}.bam
bam_sorted=${sample}_sort.bam
deduped=${sample}_deduped.bam
noMito=${sample}_noMito.bam

cores=$(more /proc/cpuinfo | grep processor | wc -l)

# move to the workind directory and begin bowtie2 alignmnet
cd ${wd}

# make the outbound directories in the folder wd
mkdir outBound || { echo "Cannot create local out bound directory. Exiting." ; exit 1 ; }
mkdir outBound/peakdeck/
mkdir outBound/peakdeck/peakdeck_75bps_10kb

#########Run bowtie2#########
bowtie2 -p ${cores} --very-sensitive --maxins 2000 -x ${bowtieindex} -1 ${fq1} -2 ${fq2} -S ${sam}

# Convert Bowtie2 sam output to bam
samtools view -1 -@ ${cores} -b -q 30 -o ${bam} ${sam}

# Sort bam file
sambamba sort -m ${memory} ${bam} -o ${bam_sorted}

# Get the chromosome sizes based ont he bam file
samtools view -H ${bam_sorted} | awk '$1 == "@SQ"{print substr($2,4)"\t"substr($3, 4)}' > chrom.sizes

# Mark duplicate reads and remove them. This was originally picard, but Java is slow so we I sambamba
sambamba markdup -r -t ${cores} -l 0 --hash-table-size=10000000 --overflow-list-size=10000000 --io-buffer-size=5000 ${bam_sorted} ${deduped}

########End and cleanup###########
# Note: the trailing "/" for "outBound/" is important! It has the effect of trimming "outBound" from the path of the files copied.
rsync -av outBound/ ${outDir}/

rm -r ${wd}