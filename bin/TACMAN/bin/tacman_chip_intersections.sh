#!/bin/bash
for fn in /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/output/BINS/*.bin.bed;
    do
    filename1=`basename ${fn}`
    for dirs in /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/output/chip/*;
        do
        dir_base=`basename ${dirs}`
        for mode in MODE1 MODE2 MODE3 MODE4;
            do
            echo $dirs/$mode
            chip_list=/Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/output/chip/"$dir_base"/*_"$mode".bed
            bedtools intersect -a $fn -b $chip_list -wb -wa > /Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/output/intersection/chip/"$filename1"_"$dir_base"_"$mode".bed
            done
        done 
    done