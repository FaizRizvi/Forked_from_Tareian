
#!/bin/bash

# You have five parameters to hardcode.
#
# The five params are:
# species - a species (id), which tells the script which reference file to use - NOTE: this remains set at "mm10" for now, so - really four params to hardcode.
# sample  - a sample (id), for the ATACseq run
# sampleName - a sampleName
# inDir      - the path to the input directory that has the fastq files
# outDir     - the path to the output directory
# 
# The †¡sample†¢ text in param (2) is used to select the appropriate fastq
# file in the input directory.
# †¡sampleName†¢ is used in the name of the output files.
# 
# †¡sample†¢ is used to select the input
# file(s, rather than †¡sampleName†¢, in case there is more than one fastq
# file to be used as input.

# ALREADY CUSTOMIZED
scriptDir="/data/miraldiLab/seqDatabase/Scripts"
echo  "scriptDir= " $scriptDir

# ALREADY CUSTOMIZED
species='hg19'

# NEED TO CUSTOMIZE
inDir='/data/miraldiLab/LCL_repo/reference_datasets/accesome/ENCODE_DNASE_GM12878_FASTQ/ENCFF001CWQ'
outDir='/data/miraldiLab/LCL_repo/20180701_ENCODE_DNASE_pipeline_results/output/ENCFF001CWQ'

# NEED TO CUSTOMIZE
sample='ENCFF001CWQ'
sampleName='ENCFF001CWQ'

# REMEMBER: create your Input subdirectory (obviously) before invocation. (The outDir subdirectory will be created automatically.)
# REMEMBER: move your input fastq.gz file(s) into the InDir before invocation, as well as the customized *.sh and *.bat files.

echo "CHECKPOINT ZERO: working in $(pwd) on $(hostname)"
echo "Settings:"
echo "species= " $species
echo "sample= " $sample
echo "sampleName= " $sampleName
echo "inDir= " $inDir
echo "outDir= " $outDir

# %%%%%%%%%%%%%%%%%%%%%%%%%%%

function log {
    echo $(date +%F_%T) $$ ${BASHPID} $1
}

log "log output with time: START"

export PATH
echo "PATH=:" $PATH

# ALREADY CUSTOMIZED
# Note: location has been hard coded
#
MarkDuplicatesJar='/usr/local/picard/1.89/jar/MarkDuplicates.jar'
echo "MarkDuplicatesJar=:" MarkDuplicatesJar

# ALREADY CUSTOMIZED
Peakdeck='/data/miraldiLab/emiraldiMariaBin/peakdeck.pl'

# ALREADY CUSTOMIZED - for use with hg19

if [[ ${species} == "hg19" ]]
then
#    BedRef='/mnt/hdfs/emiraldi/GenomeInf/mm10/Genes/complete_mouse_mm10_4col.bed'
#	acquired from UCSC table browser using hg19 genes known
     BedRef='/data/miraldiLab/LCL_repo/reference_datasets/genome/hg19/hg19_genes.bed'

#    BowtieIndex='/mnt/xfs1/bioinfo/data/Illumina/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome'
     # The six index files for Bowtie2 for mm10 were already installed on the BMI cluster, here:
#     BowtieIndex='/database/bowtie2/index/Mus_musculus/UCSC/mm10'
# 10/12/17  from Roberto:
		#6-28-2018-TAC- the bowtie index needs to point to the folder containing the files and then to the filename after the slash- the index should be .bt2....
      BowtieIndex='/database/bowtie2/index/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome'

    # TODO: EffGenSize??? Use same values as for mm9.

    echo "INFO: hg19 uses hg19 EffGenSize"
    EffGenSize=2864785220

    KgRef='/database/bowtie2/index/Homo_sapiens/UCSC/hg19/Annotation/Genes/kgXref.txt'
    RSRef='/data/miraldiLab/LCL_repo/reference_datasets/genome/hg19/hg19_GRCh37_Refseq.bed'

    echo "INFO: hg19 now using hg19/Genes/kgXref.txt rather than hg19/kgXref.txt"

    echo "Bowtie settings for species $species:"
    echo "BedRef= " $BedRef
    echo "BowtieIndex= " $BowtieIndex
    echo "EffGenSize= " $EffGenSize
    echo "KgRef= " $KgRef
    echo "RSRef= " $RSRef
else
    Usage
    exit 1
fi

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

echo "CHECKPOINT ONE: working in $(pwd) on $(hostname)"

# Create output directory if it does not already exist.
mkdir -p ${outDir} || { echo "Cannot create \"${outDir}\". Exiting." ; exit 1 ; }
# Where was I?
wwi=$(pwd)

echo "where was I (wwi) directory= " $wwi


# Create a temp directory dedicated to the processing of this sample.

# tmpDirRoot=${QC_BFSW_TMPDIRROOT:-/tmp}
  tmpDirRoot='/data/miraldiLab/Temp'

echo  "tmpDirRoot= " $tmpDirRoot


wd=$(mktemp -d ${tmpDirRoot}/${USER}_qc_bfsw_XXXXXXXXXX)

echo "working directory wd created using mktemp= " $wd

echo "CHECKPOINT TWO: working in $(pwd) on $(hostname) (have not switched to the newly created temp dir ($wd) yet "


function cleanup {
    # Try to go back to where we started from.
    cd ${wwi}

    # 11/10/17
    # Currently, for debugging, we are NOT removing the temp directory.
    # rm -rf ${wd} || echo "Check for \"${wd}\" on $(hostname)."
}
# This runs the cleanup function when the script exits (normally or
# due to an error).
trap cleanup EXIT


# Switch to the temp directory.
cd ${wd}
echo "Working in $(pwd) on $(hostname):"


echo "CHECKPOINT THREE, after switch to wd working directory: working in $(pwd) on $(hostname)"

echo "printing out var settings:"

for v in BedRef BowtieIndex EffGenSize inDir KgRef outDir RSRef sample sampleName species
do
    eval vv=\$$v
    echo "$v ${vv}"
done

# 10/22/17 
mkdir outBound || { echo "Cannot create local out bound directory. Exiting." ; exit 1 ; }

echo "Have just created the ../outBound directory under the new temp directory in $wd"


# 10/12/17
# [tay5yu@bmiclusterp2 Scripts]$ more /proc/cpuinfo | grep processor | wc -l
# 16
# [tay5yu@bmiclusterp2 Scripts]$ 
#
# cores=$(~carriero/bin/nprocNoHT)
#
# hardcoded:
#  cores='8'
#
# dynamic, found via script:
cores=$(more /proc/cpuinfo | grep processor | wc -l)
echo "cores= " $cores


# Process all of the fastq files for this sample.
shopt -s nullglob

################

# 5/7/18

declare -a arrFiles
#for file in $inDir/*_R2_      -----THIS WAS HERE, BUT GAVE ERRORS WITH UNPAIRED. TAC 6-28-2018
for file in $inDir/*
do
    arrFiles=("${arrFiles[@]}" "$file")
done

###############
echo "ARRAY:" ${arrFiles[@]}

# 5/7/18
# in this test script, we starting by trying to set fpat to match to our R1 input file:
#     041015_SS_SI_IL17pos_GCTACGCT_L003_R1_001.fastq.gz
#
# Different ways of defining fpat, all of which resulted in an empty string for fpat here
# and fqs not being filled in.
#
# fpat="${sample}_L*_R1_*.fastq.gz"
# loose_suffix="_L*_R1_*.fastq.gz"
# loose_suffix="fastq.gz"
# loose_prefix="_L003_R1_001"
#
# fpat=$sample$loose_suffix
# fpat=`041015_SS_SI_IL17pos_GCTACGCT_L*_R1_*.fastq.gz`
# fpat="041015_SS_SI_IL17pos_GCTACGCT_L003_R1_001.fastq.gz"
# fpat="$sample$loose_prefix.$loose_suffix"
#
# fqs=( ${inDir}/${fpat} )

# 5/7/18
# There were problems with defining fpat - every time Roberto P. and I tried to usefully fill it
# with a search pattern, the string emptied out, as seen in its echo statement.
# So - we went another route. We simply assumed that all the *.fastq.gz files in the
# input subdirectory for this sample were there to be used, and did not try to select
# out a subset. We simply build the arrFiles var above in a for loop, and use it
# to then fill in Nick Carriero's fqs var below.

fqs=${arrFiles[@]}

echo "--------"
echo "inDir= " $inDir
echo "sample= " $sample
# echo "loose_suffix= " $loose_suffix
# echo "fpat= " $fpat
echo "fqs= " $fqs
echo "--------"

if [[ ${#fqs[@]} == 0 ]]
then
    echo "Couldn't fid fastqs that matched default pattern, switching to a looser pattern."
    # In the new 5/7/18 code above, we have already added the inDir to the filename.
    # Hence Nick's line below is commented out.
    #   fqs=( ${inDir}/${sample}*_R1[_.]*fastq* )
    #
    if [[ ${#fqs[@]} == 0 ]]
    then
	echo "Couldn't find fastqs. Giving up."
	exit 1
    fi
    #TODO: to accommodate more flexibility here, we have to relax the
    #R1/R2 substitution below. It is not inconceivable that '_R1'
    #appears somewhere in the sample name, which would lead to
    #unintended changes to the file name. Need to clean this up.
fi

echo "Running bowtie2 and samtools sort: ${fqs[*]}"

echo "CHECKPOINT FOUR: working in $(pwd) on $(hostname)"

paired=-1
for fq in ${fqs[@]}
do
    echo "fq = " $fq

    bfq=$(basename ${fq})
    bam=${bfq/.fastq*/_bowtie2.bam}

    echo "bfq = " $bfq
    echo "bam = " $bam

    r2=${fq/_R1/_R2}

    echo "r2 = " $r2

    if [[ ${fq} != ${r2} && -e ${r2} ]]
    then
	tPaired=2
	maybe2="-1 ${fq} -2 ${r2}"
	echo "Paired"
    else
	tPaired=1
	maybe2="-U ${fq}"
	echo "Unpaired"
    fi

    if [ ${paired} == -1 ]
    then
	paired=${tPaired}
    fi
    if [ ${paired} != ${tPaired} ]
    then
	echo "WARNING: Inconsistent fastq pairing!"
	echo "WARNING: Inconsistent fastq pairing!"
    fi

    echo "maybe2 = " $maybe2
    echo "tPaired = " $tPaired
    echo "paired = " $paired

# Bowtie2 "-p"parameter definition:
#
# -p/--threads NTHREADS
#
# Launch NTHREADS parallel search threads (default: 1). Threads will
# run on separate processors/cores and synchronize when parsing reads
# and outputting alignments. Searching for alignments is highly
# parallel, and speedup is close to linear. Increasing -p increases
# Bowtie 2's memory footprint. E.g. when aligning to a human genome
# index, increasing -p from 1 to 8 increases the memory footprint by a
# few hundred megabytes. This option is only available if bowtie is
# linked with the pthreads library (i.e. if BOWTIE_PTHREADS=0 is not
# specified at build time).

echo "CHECKPOINT FIVE: working in $(pwd) on $(hostname)"

    echo "Bowtie input: ${maybe2}"
    echo "cores= " $cores
    echo "BowtieIndex= " ${BowtieIndex}
    echo "bam= " ${bam}

    { bowtie2 -p ${cores} --very-sensitive \
	--maxins 2000 \
	-x ${BowtieIndex} \
	${maybe2} | samtools view -bS -q 30 - -o ${bam} ; } 2> ${bam/.bam/.log}

    # Now sort the bam in the background and move on to the next set
    # of fastq files. We run on nodes with at least 4GB per core, so
    # we should be to spend 2GB on the sort, althought to be frank,
    # this is a bit of a guess.
    # bowtie2 -p ${cores} --very-sensitive --maxins 2000 -x ${BowtieIndex} ${maybe2} -S test.sam | samtools view -S -b test.sam > ${bam}
    ( samtools sort -m 2000000000 ${bam} ${bam/.bam/_sorted} && rm ${bam} ) &
done
wait # for sorts to finish.
echo "Done. Checking for jobs:"
jobs

echo "CHECKPOINT SIX: working in $(pwd) on $(hostname)"


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

echo "Have just created ../outBound/Stats under the $wd temp directory"

echo "CHECKPOINT SEVEN: working in $(pwd) on $(hostname)"

echo "About to invoke (altered) Python program gcdB_NJC_EM.py the first time"

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

echo "CHECKPOINT EIGHT: working in $(pwd) on $(hostname)"

echo "Removing duplicates"
deduped=${sampleName}_deduped

javaINPUT=${merged}.bam 
javaOUTPUT=${deduped}.bam 
java_METRICS_FILE=outBound/Stats/out.duplicate_metrics.${sampleName}.txt

echo "MarkDuplicatesJar= " ${MarkDuplicatesJar}
echo "javaINPUT= " $javaINPUT
echo "javaOUTPUT= " $javaOUTPUT
echo "java_METRICS_FILE= " $java_METRICS_FILE
echo "deduped= " $deduped

# Troubleshooting -- not deleting bam file
# java -Xmx1g -jar ${MarkDuplicatesJar} INPUT=${merged}.bam OUTPUT=${deduped}.bam METRICS_FILE=outBound/out.duplicate_metrics.${sampleName}.txt REMOVE_DUPLICATES=true && rm ${merged}.bam
# java -Xmx1g -jar ${MarkDuplicatesJar} INPUT=$javaINPUT OUTPUT=$javaOUTPUT METRICS_FILE=$java_METRICS_FILE REMOVE_DUPLICATES=true

  java -Xmx1g -jar ${MarkDuplicatesJar} INPUT=$javaINPUT OUTPUT=$javaOUTPUT METRICS_FILE=$java_METRICS_FILE REMOVE_DUPLICATES=true && rm $javaINPUT

echo "java call= "
echo "java -Xmx1g -jar ${MarkDuplicatesJar} INPUT=$javaINPUT OUTPUT=$javaOUTPUT METRICS_FILE=$java_METRICS_FILE REMOVE_DUPLICATES=true && rm $javaINPUT"
echo "---"

echo "About to invoke: samtools index ${deduped}.bam"	

samtools index ${deduped}.bam

echo "About to invoke the altered Python program gcdB_NJC_EM.py the second time"	


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

echo "noMito= $noMito"

echo "bamutils call:"
echo "bamutils filter ${deduped}.bam outBound/${noMito}.bam -excluderef chrM && rm ${deduped}.bam"
echo "----"
bamutils filter ${deduped}.bam outBound/${noMito}.bam -excluderef chrM && rm ${deduped}.bam


samtools index outBound/${noMito}.bam

echo "CHECKPOINT NINE: working in $(pwd) on $(hostname)"


echo "Converting bam to wig/span"
# Convert the bam file to wig/span format. This involves an implicit
# pipeline, implemented via a FIFO pipe.
#
# Strangely, bam2wig.py tries to invoke wigToBigWig, which will break in
# this context. We override PATH temporarily so that bam2wig.py will
# fail to do so (because it will not find the executable).
OLDPATH=${PATH}

NEWPATH=$(dirname $(which bam2wig.py)):/usr/bin:/bin


echo "Changing PATH to ${NEWPATH}"
echo "--"
export PATH=${NEWPATH}


wig="${noMito}.wig"

echo "wig=$wig"
echo "--"

# bam2wig.py will see this pipe as a file to which it writes its
# output. The python script that converts wig to wig/span will see it
# as input.
mkfifo FIFO${wig}


echo "About to invoke the unaltered Python program wig2wigSpan.py"
echo "python command is: python ${scriptDir}/wig2wigSpan.py < FIFO${wig} > outBound/${wig} &"
echo "--"
python ${scriptDir}/wig2wigSpan.py < FIFO${wig} > outBound/${wig} &

wwspid=$!

echo "wwspid= $wwspid"
echo "About to invoke the altered Python program bam2wig.py"
# bam2wig.py \
#	--input-file=outBound/${noMito}.bam --chromSize=chrom.sizes --out-prefix=FIFO${wig/.wig}
#
echo "python commmand is: python bam2wig.py--input-file=outBound/${noMito}.bam --chromSize=chrom.sizes --out-prefix=FIFO${wig/.wig}"
echo "--"

python ${scriptDir}/bam2wig.py \
      --input-file=outBound/${noMito}.bam --chromSize=chrom.sizes --out-prefix=FIFO${wig/.wig}



# ${scriptDir}/clearFifo.py FIFO${wig}
#
echo "About to invoke the altered Python program clearFifo.py"
echo "python command is: python ${scriptDir}/clearFifo.py FIFO${wig}"
echo "--"

python ${scriptDir}/clearFifo.py FIFO${wig}

wait ${wwspid}
    
rm FIFO${wig}

# Restore path.
echo "Restoring PATH to ${OLDPATH}"
export PATH=${OLDPATH}


echo "CHECKPOINT TEN: working in $(pwd) on $(hostname)"
echo "---"

echo "Converting wig to bigwig"
echo "command is: wigToBigWig outBound/${wig} chrom.sizes outBound/${wig/.wig/.bw}"
echo "--"

# ALREADY CUSTOMIZED
# Note: location has been hard coded by Roberto P.
/usr/local/uscstools/1.0.0/wigToBigWig outBound/${wig} chrom.sizes outBound/${wig/.wig/.bw}

echo "Generating normalized bigwig."
echo "--"
sf=$(samtools idxstats outBound/${noMito}.bam | awk '{s = s + $3}END{printf "%.12e", 1000000./s}')
echo "Scale factor: ${sf}"

# Generate the scaled bedgraph file.
#
# 
echo "bedtools command is: bedtools genomecov -ibam outBound/${noMito}.bam -bg -scale ${sf} -g chrom.sizes -split > ${noMito}.bg"
echo "--"
# bedtools genomecov -ibam outBound/${noMito}.bam -bg -scale ${sf} -g chrome.sizes -split > ${noMito}.bg
  bedtools genomecov -ibam outBound/${noMito}.bam -bg -scale ${sf} -g chrom.sizes -split > ${noMito}.bg

# Generate the bigwig for it.
echo "wigToBigWig command is: wigToBigWig -clip ${noMito}.bg chrom.sizes outBound/${noMito}_norm.bw"
echo "--"

/usr/local/uscstools/1.0.0/wigToBigWig -clip ${noMito}.bg chrom.sizes outBound/${noMito}_norm.bw


echo "CHECKPOINT ELEVEN: working in $(pwd) on $(hostname)"
echo "---"

if [[ ${RSRef} != "None" ]]
then
    ## Get a table telling you how many ATACseq bases mapped to specific gene regions

    echo "Running read_distribution.py"

    echo "About to invoke the altered Python program read_distribution.py"
    echo "python command is: python ${scriptDir}/read_distribution.py -i outBound/${noMito}.bam -r ${RSRef} > outBound/Stats/${sampleName}_geneFeatures.txt"
    echo "--"

#    read_distribution.py -i outBound/${noMito}.bam \
#	    -r ${RSRef} \
#	    > outBound/Stats/${sampleName}_geneFeatures.txt
    python ${scriptDir}/read_distribution.py -i outBound/${noMito}.bam \
	    -r ${RSRef} \
	    > outBound/Stats/${sampleName}_geneFeatures.txt


    echo "About to invoke awk script"
    echo "awk command is: awk single-quote NR >= 5 && NR <= 15{print $0} single-quote outBound/Stats/${sampleName}_geneFeatures.txt > outBound/Stats/${sampleName}_geneFeaturesTable.txt"
    echo "--"
    
    awk 'NR >= 5 && NR <= 15{print $0}' outBound/Stats/${sampleName}_geneFeatures.txt \
	    > outBound/Stats/${sampleName}_geneFeaturesTable.txt
fi


echo "CHECKPOINT TWELVE: working in $(pwd) on $(hostname)"
echo "---"

if [[ "${nopeaks+yes}" == "yes" ]]
then
    echo "Skipping peak processing."
else
    # Run two peak callers (MACS and Peakdeck) to define peaks for each
    # dataset then we will will analyze all sample peaks using
    # find_overlapping_peaks_hpc.py and sample_peak_overlaps.m to
    # visualize in a subsequent script:
    #
    # (1) peak overlap between all samples -- very useful if there are biological replicates
    #		we can look at percentage of overlapping peaks
    # (2) see how many reads fall in peaks -- this is a rough metric of signal to noise ratio 

    # Note: MACS needs to have PeakSplitter in the PATH (see above).

    # MACS IO may not be compatible with HDFS, so do this in a local directory and copy to final output directory.
    mkdir outBound/MACS
    pushd outBound/MACS || { echo "Something wrong with MACS directory." ; exit 1 ; }
    echo "Running MACS"
    echo "---"

    macs14 -t ../${noMito}.bam \
	    -n ${sampleName} \
	    -g ${EffGenSize} \
	    -B \
	    -S \
	    --call-subpeaks
	    # -g effective genome size for mm9: http://seqanswers.com/forums/showthread.php?t=4167
    popd


    # Peakdeck requires a sam file as input, so first:	

    echo "Peakdeck requires a sam file as input, so first we do this:"
    echo "samtools command is: samtools view -h -o ${sampleName}.sam outBound/${noMito}.bam"
    echo "--"

    samtools view -h -o ${sampleName}.sam outBound/${noMito}.bam

    mkdir -p outBound/peakdeck/peakdeck_75bps_10kb
    # The awk filter removes '.5' decimals from coordinates. The filter
    # will fail if the number or ordering of output fields change.

    echo "Running peakdeck.pl"
    echo "---"

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
    echo "---"
    awk '{print $0>>$1"_pd.bed"}' outBound/peakdeck/peakdeck_75bps_10kb/${sampleName}.bed
    /usr/local/uscstools/1.0.0/bigWigToWig outBound/${wig/.wig/.bw} stdout | awk '/variableStep/{l=$0 ; ofile=substr($2, 7); next}{print l"\n"$0>>ofile"_pd.wig"}'
    for pd in *_pd.bed
    do
	python ${scriptDir}/peakMaxArray.py ${pd} ${pd/.bed/.wig} > ${pd/.bed/.max} 2>outBound/${pd/_pd.bed/_pma.log}
    done | parallel -j ${cores}


#    do
#        echo "About to invoke the unaltered Python program peakMaxArray.py"
#        echo "python ${scriptDir}/peakMaxArray.py ${pd} ${pd/.bed/.wig} > ${pd/.bed/.max} 2>outBound/${pd/_pd.bed/_pma.log}"
#        echo "---"
#    done | parallel -j ${cores}


    # Merge results in same order as they appeared in the original bed
    # file. TODO: This assumes the bed does not interleave data from
    # different chroms.
    echo "Merging results in same order as they appeared in the original bed file"
    echo "---"
    ordering=$(awk '{print $1}' outBound/peakdeck/peakdeck_75bps_10kb/${sampleName}.bed | uniq)
    for c in ${ordering}
    do
	cat ${c}_pd.max
    done > outBound/peakdeck/peakdeck_75bps_10kb/${sampleName}_max.bed
fi



echo "CHECKPOINT THIRTEEN: working in $(pwd) on $(hostname)"

echo "Rsyncing results"		
# Note: the trailing "/" for "outBound/" is important! It has the
# effect of trimming "outBound" from the path of the files copied.
rsync -av outBound/ ${outDir}/

echo "CHECKPOINT FOURTEEN: DONE"

echo "Done"
log "log output with time: DONE"
