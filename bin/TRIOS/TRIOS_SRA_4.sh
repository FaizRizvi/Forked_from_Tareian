#!/bin/bash

### PARAMETERS ###

# You must have a metadata file that contains the information for how the samples are related
sample_names="/data/miraldiLab/Tareian_workspace/RNAseq_workflow/data/TRIOS_RUN_NAMES_4.txt"
index_loc="/data/miraldiLab/LCL_repo/reference_datasets/transcriptome/hg19_index"

### PIPELINE ###
for i in `cat $sample_names` 

do 
  
  mkdir $i
  cd $i  
  
  wget https://sra-download.ncbi.nlm.nih.gov/traces/sra44/SRR/004997/"$i"
  fasterq-dump ./$i

  salmon quant -i $index_loc -r "$i".fastq -o quant_"$i" --gcBias --seqBias -l A --numBootstraps 100

  rm *.fastq
  rm $i
  cd ..
done
