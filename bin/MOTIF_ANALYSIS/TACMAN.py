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
# Set up the argument parser with all of the options and defaults set up.
parser = argparse.ArgumentParser()

parser.add_argument("-b", "--BED", dest='DHS_BED', help="Input DHS Bed file", default = [], required=True)
parser.add_argument('-T','--tf', action='store', dest='tf_name_file', help='A tab-delimited file that has TF Name and Motif Name)', required=True)
parser.add_argument("-d", "--domain_bed", dest='domain_bed', help="BED file containing domain calls", required=True)
parser.add_argument("-c", "--CHIP", nargs='+', action='store', dest='CHIP_bed', help='ChIP BED files', default = [], required=True)
parser.add_argument("-G", "--META", dest='META', help="META meta file from MARIO", required=True)
parser.add_argument("-i", "--MOODS_HITS", nargs='+', action='store', dest='moods_hits', help="Motif hits from MOODS", default = [], required=True)

# set the arguments from the command line to variables in the args object
args = parser.parse_args()

# set variables from arguments
MOODS_HITS = args.moods_hits
DHS_BED = args.DHS_BED
TF_NAME_FILE = args.tf_name_file
DOMAIN_BED = args.domain_bed
CHIP_BED = args.CHIP_bed
META = args.META

##########################-----------Parameters------------##############################
# Set up styles and colors
tacman_color = colored.fg(226) + colored.attr(1)
checkpoint = stylize("################################################################################################################", tacman_color)

print (checkpoint)

# set up some reference directory file locatios
cwd = os.getcwd()
ifd = cwd + "/input"
ofd = cwd + "/output"
uwd = ofd + "/unmerged"

if not os.path.exists(ofd):
        os.makedirs(ofd)

print (stylize("TACMAN: The common working directory is: " + cwd, tacman_color))
print (stylize("TACMAN: The input directory is: " + ifd, tacman_color))
print (stylize("TACMAN: The output directory is: " + ofd, tacman_color))
print(checkpoint)

# Set up the lists that will be used
percentiles = [1, 2, 4, 10, 20, 100]
MOODS_OBJ = []
BINS = []

####################################------------CHECKPOINT 1-----------------###################################
# From here on the script will process the motif file that was produced from MOODS or provided as the -o option.
# The TF file must look ilke the following with a header column that matches the following column IDs.
# Motif_ID		TF_Name
# M00116_1.97d	TFAP2B

print (stylize("TACMAN: Parsing MOODS", tacman_color))
print(checkpoint)

os.chdir(ofd)
MOODS_wd = snippets.make_set_dir("MOODS", True)

for i in MOODS_HITS:
    length = len(MOODS_HITS)
    obj = snippets.MOODS(i, TF_NAME_FILE, ofd, DOMAIN_BED)
    obj.dict_TF_df()
    obj.parse_moods()
    obj.bed_intersect()
    obj.domain_parse()
    MOODS_OBJ.append(obj)

os.chdir(ofd)

####################################------------CHECKPOINT 2-----------------###################################
# The next step is to import the ChIP files from MARIO and arrange them by run mode, collect them into groups of common TF_start
# then take the nth percentile and return files either merged or unmerged for each tfself.
# This will also include parsing the META metadatafile from MARIOS and using it to reference.

print(checkpoint)
print (stylize("TACMAN: Parsing MARIO Data", tacman_color))

META = snippets.META(META, CHIP_BED)
META.META_parse()

chip_dir = snippets.make_set_dir("chip", True)

for i in percentiles:
    percentile = (100 - (100 / i))
    per_wd = snippets.make_set_dir("percentile_" + str(percentile), True)

    obj = snippets.MARIO(CHIP_BED, i)

    obj.parse_singles_percentile(META.unique_tf_single_df, i)

    unmerg_dir = snippets.make_set_dir("unmerged", True)

    obj.parse_replicate_percentile(META.unique_tf_rep_df, i)

    obj.merge_replicate_BEDS(META.unique_tf_rep_df, per_wd, META.unique_MODES, META.unique_TF_reps_names)

    shutil.rmtree(unmerg_dir)

    os.chdir(chip_dir)

####################################------------CHECKPOINT 3-----------------###################################
# The next step is to take the merged BED files and then create a composite BED file to BIN based on MODE

print(checkpoint)
print (stylize("TACMAN: Building reference files", tacman_color))

os.chdir(ofd)

bin_dir = snippets.make_set_dir("BINS", True)

os.chdir(ofd)

BINS = snippets.BINS(ofd, MOODS_HITS, bin_dir)

BINS.bin_group_collect(ofd, False)

BINS.bin_group_collect(ofd, True)

print (stylize("TACMAN: Intersecting binned files with ChIP files and MOODs predictions", tacman_color))

BINS.intersect_bin(META.unique_MODES)
