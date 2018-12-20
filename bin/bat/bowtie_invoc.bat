#BSUB -W 24:00
#BSUB -M 100000
#BSUB -R span[hosts=1]
#BSUB -n 16
#BSUB -J CD4_ATAC3
#BSUB -o /data/miraldiLab/Tareian/results/%J.out
#BSUB -e /data/miraldiLab/Tareian/results/%J.err

module load bowtie2/2.3.3
module load mono/git
module load peaksplitter
module load samtools/1.8.0
module load sambamba/0.6.8
module load gnuparallel
module load python/2.7.5
module load xz/5.2.3
module load openblas/0.2.20
module load jdk/1.8.0_40
module load ngsutils/0.5.9
module load gzip

memory="85G"

for i in `cat SRA_list.txt`;
do
sh bowtie2.sh $i $memory
done