#!/bin/bash

# The five params are:
# species - a species (id), which tells the script which reference file to use - NOTE: this remains set at "mm10" for now, so - really four params to hardcode.
# sample  - a sample (id), for the ATACseq run
# sampleName - a sampleName
# inDir      - the path to the input directory that has the fastq files
# outDir     - the path to the output directory
# The "sample" text in param (2) is used to select the appropriate fastq
# file in the input directory.
# "sampleName" is used in the name of the output files.

# Hardcoded Directory and file locations: Most variable names describe the item being pointed to
Peakdeck='/data/miraldiLab/emiraldiMariaBin/peakdeck.pl'
MarkDuplicatesJar='/usr/local/picard/1.89/jar/MarkDuplicates.jar'
scriptDir="/data/miraldiLab/seqDatabase/Scripts"

# Genome information
species='hg19'
BedRef='/data/miraldiLab/reference_data/genome/hg19/hg19_genes.bed'
BowtieIndex='/database/bowtie2/index/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'
EffGenSize=2864785220
KgRef='/database/bowtie2/index/Homo_sapiens/UCSC/hg19/Annotation/Genes/kgXref.txt'
RSRef='/data/miraldiLab/reference_data/genome/hg19/hg19_GRCh37_Refseq.bed'

#Input and Output directory locations
inDir='/data/miraldiLab/reference_data/Buenrostro/fastq/SRR591268'
outDir='/data/miraldiLab/reference_data/Buenrostro/fastq/SRR591268/output'
tmpDirRoot='/data/miraldiLab/Temp'

export PATH

# Sample names and name to use for folders(sampleName)
sample='SRR891268'
sampleName='SRR891268'

fq="/data/miraldiLab/reference_data/Buenrostro/fastq/SRR591268/SRR891268_1.fastq.gz"
r2="/data/miraldiLab/reference_data/Buenrostro/fastq/SRR591268/SRR891268_2.fastq.gz"

function log {
    echo $(date +%F_%T) $$ ${BASHPID} $1
}

log "log output with time: START"

echo "CHECKPOINT ZERO: working in $(pwd)"
echo "species= " $species
echo "sample= " $sample
echo "inDir= " $inDir
echo "outDir= " $outDir
echo "Bowtie settings for species $species:"
echo "BowtieIndex= " $BowtieIndex
echo "EffGenSize= " $EffGenSize

# for the bowtie index test to work, usually globbing rules have to apply, so set nullglob later in the script.
giveUp=0
for f in ${BedRef} ${BowtieIndex}* ${KgRef} ${RSRef}
do
    if [[ $f != "None" && ! -e ${f} ]]
    then
	echo "\"$f\" does not exist, giving up."
	giveUp=1
    fi
done
[[ ${giveUp} == 1 ]] && exit 1

echo "CHECKPOINT ONE: working in $(pwd)"

# Create output directory if it does not already exist and a working temp directory
mkdir -p ${outDir} || { echo "Cannot create \"${outDir}\". Exiting." ; exit 1 ; }
wwi=$(pwd)
wd=$(mktemp -d ${tmpDirRoot}/${USER}_qc_bfsw_XXXXXXXXXX)

cd ${wd}

echo "CHECKPOINT TWO: working in $(pwd)"

mkdir outBound || { echo "Cannot create local out bound directory. Exiting." ; exit 1 ; }

cores=$(more /proc/cpuinfo | grep processor | wc -l)

# Process all of the fastq files for this sample.
shopt -s nullglob

declare -a arrFiles
for file in $inDir/*_1.fastq.gz
do
    arrFiles=("${arrFiles[@]}" "$file")
done

echo "ARRAY:" ${arrFiles[@]}

echo "CHECKPOINT THREE: working in $(pwd)"

bfq=$(basename ${fq})
bam=${bfq/.fastq*/_bowtie2.bam}
bowtie_input="-1 ${fq} -2 ${r2}"

echo "fq = " $fq
echo "bam = " $bam
echo "r2 = " $r2
echo "Bowtie input: ${bowtie_input}"
echo "cores= " $cores
echo "BowtieIndex= " ${BowtieIndex}

{ bowtie2 -p ${cores} --very-sensitive --maxins 2000 -x ${BowtieIndex} \
${bowtie_input} | samtools view -bS -q 30 - -o ${bam} ; } 2> ${bam/.bam/.log}

( samtools sort -m 2000000000 ${bam} ${bam/.bam/_sorted} && rm ${bam} ) &

wait # for sorts to finish.
jobs

echo "CHECKPOINT FOUR: working in $(pwd)"
merged="${sampleName}_merged"
if [[ ${#fqs[@]} > 1 ]]
then
    echo "Running samtools merge"
    samtools merge ${merged}.bam *${sample}*_sorted.bam && rm *${sample}*_sorted.bam
else
    echo "Moving sort to merge"
    mv *${sample}*_sorted.bam ${merged}.bam
fi

echo "Running samtools index"
samtools index ${merged}.bam

# extract chrom size info from bam file.
echo "Running samtools view"
samtools view -H ${merged}.bam | awk '$1 == "@SQ"{print substr($2,4)"\t"substr($3, 4)}' > chrom.sizes

mkdir outBound/Stats

echo "CHECKPOINT FIVE: working in $(pwd)"

if [[ ${KgRef} != "None" ]]
then
    echo "Running gcdB_NJC_EM.py merged"
    python ${scriptDir}/gcdB_NJC_EM.py \
	${merged}.bam \
	${BedRef} \
	${KgRef}\
	10000 \
	outBound/Stats/${merged}_chromStats.txt

    awk 'NF == 0{exit}{print $0}' outBound/Stats/${merged}_chromStats.txt > outBound/Stats/${merged}_chromStatsTable.txt
fi

echo "CHECKPOINT SIX: working in $(pwd)"

deduped=${sampleName}_deduped
javaINPUT=${merged}.bam 
javaOUTPUT=${deduped}.bam 
java_METRICS_FILE=outBound/Stats/out.duplicate_metrics.${sampleName}.txt

echo "javaINPUT= " $javaINPUT
echo "javaOUTPUT= " $javaOUTPUT
echo "java_METRICS_FILE= " $java_METRICS_FILE
echo "deduped= " $deduped

echo "java -Xmx1g -jar ${MarkDuplicatesJar} INPUT=$javaINPUT OUTPUT=$javaOUTPUT METRICS_FILE=$java_METRICS_FILE REMOVE_DUPLICATES=true && rm $javaINPUT"
java -Xmx1g -jar ${MarkDuplicatesJar} INPUT=$javaINPUT OUTPUT=$javaOUTPUT METRICS_FILE=$java_METRICS_FILE REMOVE_DUPLICATES=true && rm $javaINPUT

echo "About to invoke: samtools index ${deduped}.bam"	
samtools index ${deduped}.bam

if [[ ${KgRef} != "None" ]]
then
    echo "Running gcdB_NJC_EM.py deduped"
	python ${scriptDir}/gcdB_NJC_EM.py \
	${deduped}.bam \
	${BedRef} \
	${KgRef} \
	10000 \
        outBound/Stats/${deduped}_chromStats.txt
    awk 'NF == 0{exit}{print $0}' outBound/Stats/${deduped}_chromStats.txt > outBound/Stats/${deduped}_chromStatsTable.txt 
fi

echo "Removing mitochondrial data"
noMito="${sampleName}_noMito"

echo "bamutils filter ${deduped}.bam outBound/${noMito}.bam -excluderef chrM && rm ${deduped}.bam"
bamutils filter ${deduped}.bam outBound/${noMito}.bam -excluderef chrM 

samtools index outBound/${noMito}.bam

echo "CHECKPOINT SEVEN: working in $(pwd) on $(hostname)"
echo "Converting bam to wig/span"
# Convert the bam file to wig/span format. This involves an implicit pipeline, implemented via a FIFO pipe.
# Strangely, bam2wig.py tries to invoke wigToBigWig, which will break in this context. 
# We override PATH temporarily so that bam2wig.py will fail to do so (because it will not find the executable).
OLDPATH=${PATH}
NEWPATH=$(dirname $(which bam2wig.py)):/usr/bin:/bin

echo "Changing PATH to ${NEWPATH}"
export PATH=${NEWPATH}

wig="${noMito}.wig"
echo "wig=$wig"

# bam2wig.py will see this pipe as a file to which it writes its output. 
# The python script that converts wig to wig/span will see it as input.
mkfifo FIFO${wig}

echo "python command is: python ${scriptDir}/wig2wigSpan.py < FIFO${wig} > outBound/${wig} &"
python ${scriptDir}/wig2wigSpan.py < FIFO${wig} > outBound/${wig} &

wwspid=$!

echo "wwspid= $wwspid"
echo "python commmand is: python bam2wig.py--input-file=outBound/${noMito}.bam --chromSize=chrom.sizes --out-prefix=FIFO${wig/.wig}"
python ${scriptDir}/bam2wig.py \
      --input-file=outBound/${noMito}.bam --chromSize=chrom.sizes --out-prefix=FIFO${wig/.wig}

echo "python command is: python ${scriptDir}/clearFifo.py FIFO${wig}"
python ${scriptDir}/clearFifo.py FIFO${wig}

wait ${wwspid}
rm FIFO${wig}

# Restore path.
echo "Restoring PATH to ${OLDPATH}"
export PATH=${OLDPATH}

echo "CHECKPOINT EIGHT: working in $(pwd)"
echo "Converting wig to bigwig"
echo "command is: wigToBigWig outBound/${wig} chrom.sizes outBound/${wig/.wig/.bw}"
/usr/local/uscstools/1.0.0/wigToBigWig outBound/${wig} chrom.sizes outBound/${wig/.wig/.bw}

echo "Generating normalized bigwig."
sf=$(samtools idxstats outBound/${noMito}.bam | awk '{s = s + $3}END{printf "%.12e", 1000000./s}')

echo "bedtools command is: bedtools genomecov -ibam outBound/${noMito}.bam -bg -scale ${sf} -g chrom.sizes -split > ${noMito}.bg"
bedtools genomecov -ibam outBound/${noMito}.bam -bg -scale ${sf} -g chrom.sizes -split > ${noMito}.bg

echo "wigToBigWig command is: wigToBigWig -clip ${noMito}.bg chrom.sizes outBound/${noMito}_norm.bw"
/usr/local/uscstools/1.0.0/wigToBigWig -clip ${noMito}.bg chrom.sizes outBound/${noMito}_norm.bw

echo "CHECKPOINT NINE: working in $(pwd)"

if [[ ${RSRef} != "None" ]]
then
    ## Get a table telling you how many ATACseq bases mapped to specific gene regions
    echo "Running read_distribution.py"
    echo "python command is: python ${scriptDir}/read_distribution.py -i outBound/${noMito}.bam -r ${RSRef} > outBound/Stats/${sampleName}_geneFeatures.txt"
    python ${scriptDir}/read_distribution.py -i outBound/${noMito}.bam \
	    -r ${RSRef} \
	    > outBound/Stats/${sampleName}_geneFeatures.txt

    echo "awk command is: awk single-quote NR >= 5 && NR <= 15{print $0} single-quote outBound/Stats/${sampleName}_geneFeatures.txt > outBound/Stats/${sampleName}_geneFeaturesTable.txt"
    awk 'NR >= 5 && NR <= 15{print $0}' outBound/Stats/${sampleName}_geneFeatures.txt \
	    > outBound/Stats/${sampleName}_geneFeaturesTable.txt
fi

echo "CHECKPOINT TEN: working in $(pwd)"

if [[ "${nopeaks+yes}" == "yes" ]]
then
    echo "Skipping peak processing."
else
    # Peakdeck requires a sam file as input, so first:	
    echo "samtools command is: samtools view -h -o ${sampleName}.sam outBound/${noMito}.bam"

    samtools view -h -o ${sampleName}.sam outBound/${noMito}.bam

    mkdir -p outBound/peakdeck/peakdeck_75bps_10kb
 
    echo "Running peakdeck.pl"
    perl ${Peakdeck} \
	    -P ${sampleName}.sam \
	    -g chrom.sizes \
	    -PVAL ON \
	    -bin 75 \
	    -STEP 25 \
	    -back 10000 \
	    -npBack 100000 \
	    | awk 'NF == 6{print $1"\t"gensub(/\.5$/, "", 1, $2)"\t"gensub(/\.5$/, "", 1, $3)"\t"$4"\t"$5"\t"$6}' > outBound/peakdeck/peakdeck_75bps_10kb/${sampleName}.bed

    echo "Find locations of maxima in peak regions"
    awk '{print $0>>$1"_pd.bed"}' outBound/peakdeck/peakdeck_75bps_10kb/${sampleName}.bed
    /usr/local/uscstools/1.0.0/bigWigToWig outBound/${wig/.wig/.bw} stdout | awk '/variableStep/{l=$0 ; ofile=substr($2, 7); next}{print l"\n"$0>>ofile"_pd.wig"}'
    for pd in *_pd.bed
    do
	python ${scriptDir}/peakMaxArray.py ${pd} ${pd/.bed/.wig} > ${pd/.bed/.max} 2>outBound/${pd/_pd.bed/_pma.log}
    done | parallel -j ${cores}

    # Merge results in same order as they appeared in the original bed
    echo "Merging results in same order as they appeared in the original bed file"
    ordering=$(awk '{print $1}' outBound/peakdeck/peakdeck_75bps_10kb/${sampleName}.bed | uniq)
    for c in ${ordering}
    do
	cat ${c}_pd.max
    done > outBound/peakdeck/peakdeck_75bps_10kb/${sampleName}_max.bed
fi

echo "CHECKPOINT ELEVEN: working in $(pwd)"

echo "Rsyncing results"		
# Note: the trailing "/" for "outBound/" is important! It has the effect of trimming "outBound" from the path of the files copied.
rsync -av outBound/ ${outDir}/

echo "Done"
log "log output with time: DONE"