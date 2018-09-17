#BSUB -W 8:00
#BSUB -M 40000
#BSUB -n 16
#BSUB -R span[hosts=1]
#BSUB -J 20180626_LCL_ENCODE_DNASE_script_run
#BSUB -o /data/miraldiLab/Tareian_workspace/ATACSeqRunStdOut/%J.out
#BSUB -e /data/miraldiLab/Tareian_workspace/ATACSeqRunStdOut/%J.err

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

sh 20180724_LCL_ENCODE_DNASE_script_hg19masked_ENCFF000SLL.sh