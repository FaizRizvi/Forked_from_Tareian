#!/bin/bash

#cd /Users/caz3so/workspaces/MOODS-python-1.9.3/scripts

#python moods_dna.py -m ./MOODS_PWM/*.pfm -s /Users/caz3so/Dropbox/thesis/data/20180812_LCL_Dnase_hg19/hg19_UNION_PEAKSET.fa -p 0.0001 -o /Users/caz3so/Dropbox/thesis/data/20180923_TACMAN_TEST/GM12878_DHS_MOTIFS_p4.txt --sep "|" --log-base 2 --batch --v

#python moods_dna.py -m ./MOODS_PWM/*.pfm -s /Users/caz3so/Dropbox/thesis/data/20180812_LCL_Dnase_hg19/hg19_UNION_PEAKSET.fa -p 0.00001 -o /Users/caz3so/Dropbox/thesis/data/20180923_TACMAN_TEST/GM12878_DHS_MOTIFS_p5.txt --sep "|" --log-base 2 --batch --v

#python moods_dna.py -m ./MOODS_PWM/*.pfm -s /Users/caz3so/Dropbox/thesis/data/20180812_LCL_Dnase_hg19/hg19_UNION_PEAKSET.fa -p 0.000001 -o /Users/caz3so/Dropbox/thesis/data/20180923_TACMAN_TEST/GM12878_DHS_MOTIFS_p6.txt --sep "|" --log-base 2 --batch --v

cd /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS

mkdir p4
cd p4

python /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/TACMAN.py -o /Users/caz3so/Dropbox/thesis/data/20180923_TACMAN_TEST/GM12878_DHS_MOTIFS_p4.txt -T /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/TF_Information_all_motifs_plus.txt -z F -d /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/Rao_GM12878_compartments_sorted.bed -r /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/GM12878_POLLARD_RE_PREDS.bed -G /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/chipseq_GM12878.geos -c /Users/caz3so/Dropbox/thesis/reference_data/CHIP/GM12878_bed_files/GM12878_bed_files/*

cd ..

mkdir p5
cd p5

python /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/TACMAN.py -o /Users/caz3so/Dropbox/thesis/data/20180923_TACMAN_TEST/GM12878_DHS_MOTIFS_p5.txt -T /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/TF_Information_all_motifs_plus.txt -z F -d /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/Rao_GM12878_compartments_sorted.bed -r /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/GM12878_POLLARD_RE_PREDS.bed -G /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/chipseq_GM12878.geos -c /Users/caz3so/Dropbox/thesis/reference_data/CHIP/GM12878_bed_files/GM12878_bed_files/*

cd ..

mkdir p6
cd p6

python /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/TACMAN.py -o /Users/caz3so/Dropbox/thesis/data/20180923_TACMAN_TEST/GM12878_DHS_MOTIFS_p6.txt -T /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/TF_Information_all_motifs_plus.txt -z F -d /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/Rao_GM12878_compartments_sorted.bed -r /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/GM12878_POLLARD_RE_PREDS.bed -G /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/chipseq_GM12878.geos -c /Users/caz3so/Dropbox/thesis/reference_data/CHIP/GM12878_bed_files/GM12878_bed_files/*

cd ..
