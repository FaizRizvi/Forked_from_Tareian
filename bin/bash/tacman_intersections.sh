#!/bin/bash

for fn in /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/output/BINS/*.bin.bed;
    do
    filename1=`basename ${fn}`
        for fn2 in /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/output/MOODS/*;
            do
            filename2=`basename ${fn2}`
            echo "Processing sample $filename1 + $filename2"
            bedtools intersect -a $fn -b $fn2 -wb -wa > /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/output/intersection/"$filename1"_"$filename2".bed
            done 
    done