#BSUB -W 24:00
#BSUB -M 16000
#BSUB -n 4
#BSUB -J PRJNA302457_Download
#BSUB -o /data/miraldiLab/reference_data/Gaublomme/%J.out
#BSUB -e /data/miraldiLab/reference_data/Gaublomme/%J.err

module load sratoolkit/2.9.0
module load gzip

## To invoke this batch job on the BMI cluster:
## 1) log into the cluster
## 2) issue this command at the command prompt:
##    bsub < <the_name_of_this_bat_file>.bat

sh PRJNA302457_fastqdump.sh

