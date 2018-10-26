#!/usr/bin/env python
from __future__ import print_function
import time
start_time = time.time()

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

# set variables from arguments
MOODS_HITS = args.moods_hits
DHS_BED = args.DHS_BED
TF_NAME_FILE = args.tf_name_file
DOMAIN_BED = args.domain_bed
CHIP_BED = args.CHIP_bed
META = args.META
BLACKLIST = args.blacklist

##########################-----------Parameters------------##############################
"""Set up the colors to be used on the display and also what the checpoint borders would look ilke"""
tacman_color = colored.fg(226) + colored.attr(1)
checkpoint = stylize("################################################################################################################", tacman_color)
print (checkpoint)

# set up some reference directory file locatios
cwd = os.getcwd()
ifd = cwd + "/input"
ofd = cwd + "/output"
uwd = ofd + "/unmerged"

# If the path does not exists for the output directory, create it
if not os.path.exists(ofd):
        os.makedirs(ofd)

# Print the names of the directories for the user
print (stylize("TACMAN: The common working directory is: " + cwd, tacman_color))
print (stylize("TACMAN: The input directory is: " + ifd, tacman_color))
print (stylize("TACMAN: The output directory is: " + ofd, tacman_color))
print (checkpoint)

# Set up the various lists that will be used. The first list is the percentile cutoff that wille be used
# the next list will hold the MOODS objects that are created
# the last list has the BIN files that will be created. these might be replaced by objects in later versions
percentiles = [1, 2, 4, 10, 20, 100]
MOODS_OBJ = []
BINS = []

####################################------------CHECKPOINT 1-----------------###################################
"""From here on the script will process the motif file that was produced from MOODS or provided as the -o option.
The TF file must look ilke the following with a header column that matches the following column IDs"""

# Motif_ID		TF_Name
# M00116_1.97d	TFAP2B

#Find the time that has elasped between the loading of the program and decide how to display it
elapsed = (time.time() - start_time)

if elapsed <= 60:
        print (stylize("---It took TACMAN %s seconds to boot up and set up---" % elapsed, tacman_color))
        print (checkpoint)

else:
        minutes = elapsed/60
        print (stylize("---It took TACMAN %s minutes to boot up and set up---" % minutes, tacman_color))
        print (checkpoint)
     
print (stylize("TACMAN: Parsing MOODS", tacman_color))

# Change the working directory to the output file directory
os.chdir(ofd)

# Make the subdirectory MOODS in the OFD
MOODS_wd = snippets.make_set_dir("MOODS", True)

# Start the moods timer for profiling
MOODS_start_time = time.time()

# For every MOODS file with a different p-value create an object
for i in MOODS_HITS:
    length = len(MOODS_HITS)
    obj = snippets.MOODS(i, TF_NAME_FILE, ofd, DOMAIN_BED)
    obj.dict_TF_df()
    obj.parse_moods()
    obj.bed_intersect()
    obj.domain_parse()
    MOODS_OBJ.append(obj)

# calculate the time it took to run the function
snippets.clock(MOODS_start_time, start_time, "Parse MOODS")
print (checkpoint)

MARIO_start_time = time.time()

# change the working directory back to the ofd
os.chdir(ofd)

####################################------------CHECKPOINT 2-----------------###################################
"""The next step is to import the ChIP files from MARIO then take the nth percentile and return merged 
This will also include parsing the META metadatafile from MARIOS and using it to reference."""
print (stylize("TACMAN: Parsing MARIO Data", tacman_color))

# Create the META file object for the ChIP data 
META = snippets.META(META, CHIP_BED)
META.META_parse()

# Output the time
snippets.clock(MARIO_start_time, start_time, "Parse MARIO")
print (checkpoint)

# make the Chip sub directory
chip_dir = snippets.make_set_dir("chip", True)

# start the time for the parsing function
Percentile_start_time = time.time()

print (stylize("TACMAN: Parsing Percentiles from ChIP bed files", tacman_color))

percentile_folder = []

# For each percentile cutoff parse the data and save in a subfolder of ChIP
for i in percentiles:
    percentile = (100 - (100 / i))
    per_wd = snippets.make_set_dir("percentile_" + str(percentile), True)

    obj = snippets.MARIO(CHIP_BED, i, BLACKLIST)

    obj.parse_singles_percentile(META.unique_tf_single_df, i)

    unmerg_dir = snippets.make_set_dir("unmerged", True)

    obj.parse_replicate_percentile(META.unique_tf_rep_df, i)
    obj.merge_replicate_BEDS(META.unique_tf_rep_df, per_wd, META.unique_MODES, META.unique_TF_reps_names)

    shutil.rmtree(unmerg_dir)
    percentile_folder.append(per_wd)

    os.chdir(chip_dir)

#Report the time for the percentile parse
snippets.clock(Percentile_start_time, start_time, "Parse Percentiles")
print (checkpoint)

####################################------------CHECKPOINT 3-----------------###################################
"""The next step is to take the merged BED files and then create a composite BED file to BIN based on MODE
Go back to the base directory"""
os.chdir(ofd)

# Create the bin directory under the ofd
bin_dir = snippets.make_set_dir("BINS", True)

# Go back to OFD since it where the file reference from
os.chdir(ofd)

# Start the timer for the Binning function
BINS_start_time = time.time()

# Create the BIN object
BINS = snippets.BINS(ofd, DHS_BED, bin_dir)

# Bin the DHS areas only. This will also remove the blacklist areas from the DHS file. 
print (stylize("TACMAN: Binnning DHS", tacman_color))
BINS.bin_group_collect(ofd, True, BLACKLIST)

# Bin the union of all files 
print (stylize("TACMAN: Binnning Union", tacman_color))
BINS.bin_group_collect(ofd, False, BLACKLIST)

# Print out the time it took to Bin the Genomes
print (checkpoint)
snippets.clock(BINS_start_time, start_time, "Binning Genomes")
print (checkpoint)

####################################------------CHECKPOINT 4-----------------###################################
"""Intersect all of the files that have been collected with the two BIN reference files"""
os.chdir(ofd)

print (stylize("TACMAN: Intersecting binned files with ChIP files and MOODs predictions", tacman_color))

intersect_dir = snippets.make_set_dir("intersection", True)

# Intersect the CHIP files with the two bin reference files
BINS.intersect_chip_bin(META.unique_MODES, BINS.dhs, percentile_folder, intersect_dir, "DHS")
BINS.intersect_chip_bin(META.unique_MODES, BINS.union, percentile_folder, intersect_dir, "UNION")

# Intersect the MOODS files with the two bin reference files
BINS.intersect_moods_bin(BINS.dhs, MOODS_wd, intersect_dir, "DHS")
BINS.intersect_moods_bin(BINS.union, MOODS_wd, intersect_dir, "UNION")

####################################------------CHECKPOINT 5-----------------###################################
"""Start analyzing the intersection data gathered. This will include grouping the datae, 
parsing a lot, and then calculating P/R. I will also do all analysis in this section that 
I want to perform. """

