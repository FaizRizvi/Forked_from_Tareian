#BSUB -W 24:00
#BSUB -M 75000
#BSUB -n 16
#BSUB -J Buenrostro_1
#BSUB -o /data/miraldiLab/reference_data/Buenrostro/%J.out
#BSUB -e /data/miraldiLab/reference_data/Buenrostro/%J.err

module load bedtools/2.17.0
module load bowtie2/2.3.3
module load igvtools/2.4.0
module load mono/git
module load peaksplitter
module load samtools/0.1.19

module load gnuparallel
module load picard/1.89
module load python/2.7.10

module load xz/5.2.3
module load openblas/0.2.20

module load jdk/1.8.0_40
module load ngsutils/0.5.9
module load gzip

sh 20180724_GM12878_Buenrostro_SRR891268.sh
