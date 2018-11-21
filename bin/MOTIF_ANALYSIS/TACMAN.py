#!/usr/bin/env python
from __future__ import print_function

import snippets
import os
import sys
import argparse
import pandas as pd
import seaborn as sns
import numpy as np
import colored
from itertools import groupby, chain
from colored import stylize
import shutil

##########################-----------ARGUMENT PARSER------------##############################
"""Set up the argument parser with all of the options and defaults set up."""
parser = argparse.ArgumentParser()

parser.add_argument("-b", "--BED", dest='DHS_BED', help="Input DHS Bed file", required=True)
parser.add_argument('-T','--tf', action='store', dest='tf_name_file', help='A tab-delimited file that has TF Name and Motif Name)', required=True)
parser.add_argument("-d", "--domain_bed", dest='domain_bed', help="BED file containing domain calls", required=True)
parser.add_argument("-c", "--CHIP", nargs='+', action='store', dest='CHIP_bed', help='ChIP BED files', default = [], required=True)
parser.add_argument("-G", "--META", dest='META', help="META meta file from MARIO", required=True)
parser.add_argument("-i", "--MOODS_HITS", nargs='+', action='store', dest='moods_hits', help="Motif hits from MOODS", default = [], required=True)
parser.add_argument("-bl", "--blacklist", dest='blacklist', help="Blacklist regions", required=True)

# set the arguments from the command line to variables in the args object
args = parser.parse_args()

##########################-----------Parameters------------##############################
"""Set up the colors to be used on the display and also what the checpoint borders would look ilke
Also set up the directories. This needs to be changed in the future to makte it easy to set input and output paths"""
tacman_color = colored.fg(226) + colored.attr(1)

# set up some reference directory file locatios
CWD = os.getcwd()
OFD = CWD + "/output"

# If the path does not exists for the output directory, create it
if not os.path.exists(OFD):
        os.makedirs(OFD)

# Print the names of the directories for the user
print (stylize("TACMAN: The common working directory is: " + CWD, tacman_color))
print (stylize("TACMAN: The output directory is: " + OFD, tacman_color))

# Set up the various lists that will be used. The first list is the percentile cutoff that wille be used
percentiles = [1, 2, 4, 10, 20, 100]

#These are a bunch of lists that I will populate. 
# the next list will hold the MOODS objects that are created
# the last list has the BIN files that will be created. these might be replaced by objects in later versions
percentile_len = len(percentiles)
percentile_counter = 1
MOODS_OBJ = []
BINS = []
percentile_folder = []
MOODS_len = len(args.moods_hits)
MOODS_counter = 1

####################################------------CHECKPOINT 1-----------------###################################
"""From here on the script will process the motif file that was produced from MOODS or provided as the -o option.
The TF file must look ilke the following with a header column that matches the following column IDs

Motif_ID	TF_Name
M00116_1.97d	TFAP2B
 """   
print (stylize("TACMAN: Parsing MOODS", tacman_color))

# Change the working directory to the output file directory
os.chdir(OFD)

# Make the subdirectory MOODS in the output file directory
MOODS_wd = snippets.make_set_dir("MOODS", True)

# Create META object
META = snippets.META(args.META, args.CHIP_bed, args.tf_name_file)

# For every MOODS file with a different p-value create an object
for i in args.moods_hits:
    print ("working on file: " + str(MOODS_counter) + "/" + str(MOODS_len))
    obj = snippets.MOODS(i, OFD, args.domain_bed, META.tf_dict, META.unique_TFS)
    MOODS_counter = MOODS_counter + 1
    MOODS_OBJ.append(obj)

# change the working directory back to the OFD
os.chdir(OFD)

####################################------------CHECKPOINT 2-----------------###################################
"""The next step is to import the ChIP files from MARIO then take the nth percentile and return merged 
This will also include parsing the META metadatafile from MARIOS and using it to reference."""

# make the Chip sub directory
chip_dir = snippets.make_set_dir("chip", True)

print (stylize("TACMAN: Parsing Percentiles from ChIP bed files", tacman_color))

# For each percentile cutoff parse the data and save in a subfolder of ChIP
for i in percentiles:
    percentile = (100 - (100 / i))
    
    print ("Working on parsing percentile cutoff: " + str(percentile_counter) + "/" + str(percentile_len))
    per_wd = snippets.make_set_dir("percentile_" + str(percentile), True)

    obj = snippets.MARIO(args.CHIP_bed, i, args.blacklist)

    obj.parse_singles_percentile(META.unique_tf_single_df, i)
    
    unmerg_dir = snippets.make_set_dir("unmerged", True)

    obj.parse_replicate_percentile(META.unique_tf_rep_df, i, per_wd, META.unique_MODES, META.unique_TF_reps_names)

    shutil.rmtree(unmerg_dir)

    percentile_folder.append(per_wd)
    percentile_counter = percentile_counter + 1
    
    os.chdir(chip_dir)

####################################------------CHECKPOINT 3-----------------###################################
"""The next step is to take the merged BED files and then create a composite BED file to BIN based on MODE
Go back to the base directory"""
os.chdir(OFD)

# Create the bin directory under the OFD
bin_dir = snippets.make_set_dir("BINS", True)

# Go back to OFD since it where the file reference from
os.chdir(OFD)

# Create the BIN object
BINS = snippets.BINS(chip_dir, args.DHS_BED, bin_dir, OFD, chip_dir)

# Bin the DHS areas only. This will also remove the blacklist areas from the DHS file. 
print (stylize("TACMAN: Binnning DHS", tacman_color))
BINS.bin_group_collect(OFD, True, args.blacklist)

# Bin the union of all files 
print (stylize("TACMAN: Binnning Union", tacman_color))
BINS.bin_group_collect(OFD, False, args.blacklist)

####################################------------CHECKPOINT 4-----------------###################################
"""Intersect all of the files that have been collected with the two BIN reference files"""
os.chdir(OFD)

print (stylize("TACMAN: Intersecting binned files with ChIP files", tacman_color))
intersect_dir = snippets.make_set_dir("intersection", True)
intersect_chip_dir = snippets.make_set_dir("CHIP", True)

# Intersect the CHIP files with the two bin reference files
print (stylize("TACMAN: Intersecting binned DHS with ChIP files", tacman_color))
BINS.intersect_chip_bin(META.unique_MODES, BINS.dhs, percentile_folder, intersect_chip_dir, "DHS", META.unique_TF_reps_names)

print (stylize("TACMAN: Intersecting binned UNION with ChIP files", tacman_color))
BINS.intersect_chip_bin(META.unique_MODES, BINS.union, percentile_folder, intersect_chip_dir, "UNION", META.unique_TF_reps_names)

os.chdir(intersect_dir)

intersect_MOODS_dir = snippets.make_set_dir("MOODS", True)

# Intersect the MOODS files with the two bin reference files
print (stylize("TACMAN: Intersecting binned DHS with MOODs predictions", tacman_color))
BINS.intersect_moods_bin(BINS.dhs, MOODS_wd, intersect_MOODS_dir, "DHS", META.unique_TF_reps_names)

print (stylize("TACMAN: Intersecting binned UNION with MOODs predictions", tacman_color))
BINS.intersect_moods_bin(BINS.union, MOODS_wd, intersect_MOODS_dir, "UNION", META.unique_TF_reps_names)

####################################------------CHECKPOINT 5-----------------###################################
"""Start analyzing the intersection data gathered. This will include grouping the datae, 
parsing a lot, and then calculating P/R. I will also do all analysis in this section that 
I want to perform. """

print (stylize("TACMAN: Analyzing Results", tacman_color))

os.chdir(OFD)

results_dir = snippets.make_set_dir("results", True)

