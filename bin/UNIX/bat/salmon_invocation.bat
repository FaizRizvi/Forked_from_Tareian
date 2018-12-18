#BSUB -W 24:00
#BSUB -M 36000
#BSUB -n 16
#BSUB -R span[hosts=1]
#BSUB -J salmon_cluster_quant_test
#BSUB -o /data/miraldiLab/seqDatabase/ATACSeqRunStdOut/%J.out
#BSUB -e /data/miraldiLab/seqDatabase/ATACSeqRunStdOut/%J.err

module load python/2.7.10
module load gzip
module load salmon

## To invoke this batch job on the BMI cluster:
## 1) log into the cluster
## 2) issue this command at the command prompt:
##    bsub < <the_name_of_this_bat_file>.bat

sh wget_md5sum_salmon.sh

