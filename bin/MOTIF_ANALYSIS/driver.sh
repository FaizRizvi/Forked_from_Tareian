#!/bin/bash
shopt -s nullglob

PWM=/Users/caz3so/Dropbox/thesis/reference_data/MOTIF_SCAN/MOODS_PWM/*



python ./TACMAN.py -o test.txt -m $PWM \
-s /Users/caz3so/Dropbox/thesis/data/20180812_LCL_Dnase_hg19/hg19_UNION_PEAKSET.fa \
--sep "|" --batch -T TF_Information_all_motifs_plus.txt -z T -p .0001 \
-d Rao_GM12878_compartments_sorted.bed -r GM12878_POLLARD_RE_PREDS.bed -G chipseq_GM12878.geos \
-c /Users/caz3so/Dropbox/thesis/reference_data/CHIP/GM12878_bed_files/GM12878_bed_files/*
