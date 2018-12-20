#!/bin/bash

#########Input#########
sample=$1
memory=$2

outDirRoot=/data/miraldiLab/Tareian/results
tmpDirRoot=/data/miraldiLab/Temp
bowtieindex="/home/tacazares/workspace/bowtie2/hg19/hg19"

#########Set Parameters from input#########
outDir=${outDirRoot}/${sample}
mkdir -p ${outDir} || { echo "Cannot create \"${outDir}\". Exiting." ; exit 1 ; }
wwi=$(pwd)
wd=$(mktemp -d ${tmpDirRoot}/${USER}_qc_bfsw_XXXXXXXXXX)

fq1=${sample}_1.fastq
fq2=${sample}_2.fastq

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
bowtie2 -p ${cores} --very-sensitive --maxins 2000 -x /home/tacazares/workspace/bowtie2/hg19/hg19 -1 /home/tacazares/workspace/fastq/Buenrostro/SRR891268_1.fastq -2 /home/tacazares/workspace/fastq/Buenrostro/SRR891268_2.fastq -S ${sam}

# Convert Bowtie2 sam output to bam
samtools view -1 -@ ${cores} -b -q 30 -o ${bam} ${sam}

# Sort bam file
sambamba sort -m ${memory} ${bam} -o ${bam_sorted}

# Get the chromosome sizes based ont he bam file
samtools view -H ${bam_sorted} | awk '$1 == "@SQ"{print substr($2,4)"\t"substr($3, 4)}' > chrom.sizes

# Mark duplicate reads and remove them. This was originally picard, but Java is slow so we I sambamba
sambamba markdup -r -t ${cores} -l 0 --hash-table-size=10000000 --overflow-list-size=10000000 --io-buffer-size=5000 ${bam_sorted} ${deduped}

# Remove the mitochondrial reads
samtools idxstats ${deduped} | cut -f 1 | grep -v MT | xargs samtools view -b ${deduped} > outBound/${noMito}

# Index the no mito BAM
samtools index outBound/${noMito}

#########Peakdeck#############
# Peakdeck requires a sam file as input, so first:	
samtools view -@ ${cores} -h -o ${sam} outBound/${noMito}

# Run Peakdeck. The awk filter removes '.5' decimals from coordinates. The filterwill fail if the number or ordering of output fields change.
perl ${Peakdeck} -P ${sam} -g chrom.sizes -PVAL ON -bin 75 -STEP 25 -back 10000 -npBack 100000 \
    | awk 'NF == 6{print $1"\t"gensub(/\.5$/, "", 1, $2)"\t"gensub(/\.5$/, "", 1, $3)"\t"$4"\t"$5"\t"$6}' > outBound/peakdeck/peakdeck_75bps_10kb/${sample}.bed

#########MACS2#############
# The only required input for MACS2 is the treatment file. I wanted to get as close to Peakdeck as possible
MACS2 -t outBound/${noMito} -n ${sample} -o outbound/MACS2 --nomodel --nolambda -f BAM -g ${EffGenSize}

########End and cleanup###########
# Note: the trailing "/" for "outBound/" is important! It has the effect of trimming "outBound" from the path of the files copied.
rsync -av outBound/ ${outDir}/

rm -r ${wd}