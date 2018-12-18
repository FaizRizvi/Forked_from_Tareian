#BSUB -W 48:00
#BSUB -M 16000
#BSUB -n 8
#BSUB -R span[hosts=1]
#BSUB -J TRIOS_FASTQDUMP
#BSUB -o /data/miraldiLab/Tareian_workspace/ATACSeqRunStdOut/%J.out
#BSUB -e /data/miraldiLab/Tareian_workspace/ATACSeqRunStdOut/%J.err

module load salmon
module load 
module load sratoolkit/2.9.0

## To invoke this batch job on the BMI cluster:
## 1) log into the cluster
## 2) issue this command at the command prompt:
##    bsub < <the_name_of_this_bat_file>.bat

sh wget_md5sum_salmon.sh

-bash-4.1$ bsub -W 12:00 -M 32000 -n 8 -Is bash
