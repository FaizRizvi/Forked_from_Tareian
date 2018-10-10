#!/usr/bin/env python
from __future__ import print_function

import snippets
import os
import sys
import argparse
import pandas as pd
import seaborn as sns
import pybedtools
import math
import numpy as np
import scipy
import colored
from itertools import groupby, chain
from colored import stylize
import fish

bird = fish.Bird()

##########################-----------ARGUMENT PARSER------------##############################
#Set up the argument parser with all of the options and defaults set up.
################################################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--BED",  help="Input DHS Bed file", default = [], required=True)
parser.add_argument('-T','--tf', action='store', dest='tf_name_file', help='A tab-delimited file that has TF Name and Motif Name)', required=True)
parser.add_argument("-d", "--domain_bed", dest='domain_bed', help="BED file containing domain calls", required=True)
parser.add_argument("-r", "--RE_bed", dest='RE_bed', help="BED file containing RE annotations", required=True)
parser.add_argument("-c", "--CHIP", nargs='+', action='store', dest='CHIP_bed', help='ChIP BED files', default = [], required=True)
parser.add_argument("-G", "--GEOS_meta", dest='GEOS_meta', help="GEOS meta file from MARIO", required=True)
parser.add_argument("-i", "--MOODS_HITS", nargs='+', action='store', dest='moods_hits', help="Motif hits from MOODS", default = [], required=True)
parser.add_argument("-MP", "--MP_bool", dest='mp_bool', help="T or F for Moods Parse", required=True)

#set the arguments from the command line to variables in the args object
args = parser.parse_args()

#set variables from arguments
MOODS_HITS = args.moods_hits
TF_NAME_FILE = args.tf_name_file
DOMAIN_BED = args.domain_bed
RE_BED = args.RE_bed
CHIP_BED = args.CHIP_bed
GEOS_META = args.GEOS_meta
MP_BOOL = args.mp_bool

#Styles
tacman_color = colored.fg(226) + colored.attr(1)

checkpoint = stylize("################################################################################################################", tacman_color)
print (checkpoint)

#set up some reference directory file locatios
cwd = os.getcwd()
ifd = cwd + "/input"
ofd = cwd + "/output"

print (stylize("TACMAN: The common working directory is: " + cwd, tacman_color))
print (stylize("TACMAN: The input directory is: " + ifd, tacman_color))

if not os.path.exists(ofd):
        os.makedirs(ofd)

print (stylize("TACMAN: The output directory is: " + ofd, tacman_color))
print(checkpoint)

####################################------------CHECKPOINT 1-----------------###################################
# From here on the script will process the motif file that was produced from MOODS or provided as the -o option.
################################################################################################################
# The TF file must look something ilke the following with a header column that matches the following column IDs.
#
# Motif_ID		TF_Name
# M00116_1.97d	TFAP2B

MOODS_BED_LIST = []

if MP_BOOL == "T":
    print (stylize("TACMAN: Parsing MOODS", tacman_color))
    print(checkpoint)
    for i in MOODS_HITS:
        name = snippets.parse_moods(i, TF_NAME_FILE, ofd)
        MOODS_BED_LIST.append(name)
else:
	print(stylize("TACMAN: Not Parsing MOODS", tacman_color))

####################################------------CHECKPOINT 2-----------------###################################
# From here the parsed Moods file will be in a BED format and will be intersected with the RE lists of interest and also the domains.
################################################################################################################
print(checkpoint)

MOTIF_PATH = []

for j in MOODS_BED_LIST:
    print (stylize("TACMAN: Intersecting Domain Bed With: " + j, tacman_color))
    print(checkpoint)

    TF_domain_intersect = snippets.bed_intersect(j, DOMAIN_BED, True)

    print(checkpoint)
    print (stylize("TACMAN: Parsing Domain BED intersected with: " + j, tacman_color))
    print(checkpoint)

    MOODS_HITS_BED_PATHS = snippets.domain_parse(TF_domain_intersect, j)

    MOTIF_PATH = MOTIF_PATH.append(MOODS_HITS_BED_PATHS)

####################################------------CHECKPOINT 4-----------------###################################
# The next step is to import the ChIP files from MARIO and arrange them by run mode, collect them into groups of common TF_start
# then take the nth percentile and return files either merged or unmerged for each tfself.
# This will also include parsing the GEOS metadatafile from MARIOS and using it to reference.
################################################################################################################
print(checkpoint)
print (stylize("TACMAN: Parsing MARIO Data", tacman_color))
print(checkpoint)

CHIP_df = snippets.parse_MARIO(GEOS_META, CHIP_BED)

print(checkpoint)
print (stylize("TACMAN: Making ChIP folder", tacman_color))
print(checkpoint)

snippets.make_set_dir("chip")

print(checkpoint)
print (stylize("TACMAN: Parsing out the 75th percentile", tacman_color))
print(checkpoint)

CHIP_df["file_path"].apply(snippets.parse_percentile)

print(checkpoint)
print (stylize("TACMAN: Making BED folder", tacman_color))
print(checkpoint)

snippets.make_set_dir("merged")

print(checkpoint)
print (stylize("TACMAN: Merging replicate TF BEDS", tacman_color))
print(checkpoint)

snippets.parse_replicate_BEDS(CHIP_df, MOTIF_PATH)

####################################------------CHECKPOINT 5-----------------###################################
# The next step is to take the merged BED files and then create a composite BED file to BIN based on MODE
################################################################################################################
print(checkpoint)
print (stylize("TACMAN: Binning genome based on MODE", tacman_color))
print(checkpoint)

print (stylize("TACMAN: Building reference files", tacman_color))

MODES = ["MODE1", "MODE2", "MODE3", "MODE4"]
MODES = pd.DataFrame(MODES)

#In this line of code you will break down the different MODES into sets that will be binned. The next step is sorting.
MODES[0].apply(snippets.bin_group_collect, args = (CHIP_group_df, MOTIF_PATH, files_list))

####################################------------CHECKPOINT 6-----------------###################################
# The next step is to take the merged BED files and then create a composite BED file to BIN based on MODE
################################################################################################################
print (stylize("TACMAN: Intersecting binned files with ChIP files and MOODs predictions", tacman_color))

snippets.intersect_bin(MOTIF_PATH)
