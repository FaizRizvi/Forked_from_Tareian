#BSUB -W 12:00
#BSUB -M 250000
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -J sleuth_cluster_quant_test
#BSUB -o /data/miraldiLab/Tareian_workspace/DEseq/%J.out
#BSUB -e /data/miraldiLab/Tareian_workspace/DEseq/%J.err

module load R/3.5.0

## To invoke this batch job on the BMI cluster:
## 1) log into the cluster
## 2) issue this command at the command prompt:
##    bsub < <the_name_of_this_bat_file>.bat

R -f sleuth_wasabi.r

