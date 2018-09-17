#BSUB -W 4:00
#BSUB -M 100000
#BSUB -n 16
#BSUB -R span[hosts=1]
#BSUB -J salmon_cluster_quant_test
#BSUB -o /data/miraldiLab/seqDatabase/ATACSeqRunStdOut/%J.out
#BSUB -e /data/miraldiLab/seqDatabase/ATACSeqRunStdOut/%J.err

### this script runs through how to run a live session on the cluster
### and perform DEseq Analysis

module load R/3.5.0

R -f DESEQ2_LAB.r