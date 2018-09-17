#!/bin/bash

### PARAMETERS ###
# The parameters that are needed are:

#1: an indexed transcriptome
index_loc="/data/miraldiLab/LCL_repo/reference_datasets/hg19_index"

#2: a file of md5sums
md5sum_file_loc="/data/miraldiLab/Tareian_workspace/RNAseq_workflow/geuvadis_md5_fastq.m5"

#3: list of the sample names
sampleFile="/data/miraldiLab/Tareian_workspace/RNAseq_workflow/sample_names_unique_batch4.txt"

### PIPELINE ###
for i in `cat $sampleFile`
do
  mkdir $i
  cd $i

  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/"$i"/"$i"_1.fastq.gz
  wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/"$i"/"$i"_2.fastq.gz

  # GREP the names from the metafile and then create a temp file to use md5sum
  grep "$i"_1.fastq.gz $md5sum_file_loc > temp_hash.m5
  grep "$i"_2.fastq.gz $md5sum_file_loc >> temp_hash.m5

  md5sum -c temp_hash.m5

  salmon quant -i $index_loc -1 "$i"_1.fastq.gz -2  "$i"_2.fastq.gz -o quant_"$i" --gcBias --seqBias -l A --numBootstraps 100

  rm "$i"_1.fastq.gz
  rm "$i"_2.fastq.gz
  rm temp_hash.m5
  cd ..
done
