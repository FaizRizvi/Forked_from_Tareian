#!/usr/bin/env python
from __future__ import print_function

import MOODS.scan
import MOODS.tools
import MOODS.parsers
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
import time
from tqdm import tqdm

from itertools import groupby, chain

##########################-----------ARGUMENT PARSER------------##############################
#Set up the argument parser with all of the options and defaults set up.
################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbosity", action="count", default=0, help='verbose (-vv, -vvv for more)')
parser.add_argument("-b", "--BED",  help="Input DHS Bed file")
parser.add_argument("-z", "--MOODS", help="Run Moods", action="store")
parser.add_argument('-T','--tf', action='store', dest='tf_name_file', help='A tab-delimited file that has TF Name and Motif Name)')
parser.add_argument("-d", "--domain_bed", dest='domain_bed', help="BED file containing domain calls")
parser.add_argument("-r", "--RE_bed", dest='RE_bed', help="BED file containing RE annotations")
parser.add_argument("-c", "--CHIP", nargs='+', action='store', dest='CHIP_bed', help='ChIP BED files', default = [])
parser.add_argument("-G", "--GEOS_meta", dest='GEOS_meta', help="GEOS meta file from MARIO")

#MOODS OPTIONS - input files
input_group = parser.add_argument_group("input files (at least one matrix and sequence file required)")
input_group.add_argument('-m','--matrices', metavar='M', nargs='+', action='store', dest='matrix_files', help='matrix files (count/frequency, will be converted to PWM before matching)', default = [])
input_group.add_argument('-S','--score-matrices', metavar='M', nargs='+', action='store', dest='lo_matrix_files', help='matrix files (PWM/other score matrix, will be matched directly)', default = [])
input_group.add_argument('-s','--sequences', metavar='S', nargs='+', action='store', dest='sequence_files', help='sequence files', default = [])

#MOODS OPTIONS - thresholds
th_group = parser.add_argument_group("threshold selection (exactly one required)")
th_group.add_argument('-p','--p-value', metavar='p', action='store', dest='p_val', type=float, help='compute threshold from p-value')
th_group.add_argument('-t','--threshold', metavar='T', action='store', dest='t', type=float, help='use specified absolute threshold')
th_group.add_argument('-B','--best-hits', metavar='n', action='store', dest='max_hits', type=int, help='return at least the specified amount of best matches')

#MOODS OPTIONS - output
out_group = parser.add_argument_group("output (optional)")
out_group.add_argument('-o','--outfile', metavar='outfile', action='store', dest='output_file', help='output to file instead of standard output', required=True)
out_group.add_argument('--sep', metavar='S', action='store', dest='sep', default=",", help='set field separator in out (default ",")')

#MOODS OPTIONS - behaviour
option_group = parser.add_argument_group("search and model behaviour (optional)")
option_group.add_argument('-R','--no-rc', dest='no_rc', action="store_true", help='disable matching versus reverse complement strand')
option_group.add_argument('--no-snps', dest='no_snps', action="store_true", help='ignore IUPAC symbols coding multiple nucleotides')
option_group.add_argument('--batch', dest='batch', action="store_true", help='do not recompute thresholds from p-values for each sequence separately (recommended when dealing with lots of short sequences)', default = "batch")
option_group.add_argument('--bg', metavar=('pA', 'pC', 'pG', 'pT'), nargs=4, action='store', type=float, dest='bg', default = [0.25,0.25,0.25,0.25], help='background distribution for computing thresholds from p-value with --batch (default is 0.25 for all alleles)')
option_group.add_argument('--ps', metavar='p', action='store', dest='ps', type=float, help='specify pseudocount for log-odds conversion (default = 0.1)', default = 0.01)# bg
option_group.add_argument('--log-base', metavar='x', action='store', dest='log_base', type=float, help='logarithm base for log-odds conversion (default natural logarithm)')
option_group.add_argument('--lo-bg', metavar=('pA', 'pC', 'pG', 'pT'), nargs=4, action='store', type=float, dest='lo_bg', default = [0.25,0.25,0.25,0.25], help='background distribution for log-odds conversion (default is 0.25 for all alleles)')

#set the arguments from the command line to variables in the args object
args = parser.parse_args()

#set variables from arguments
MOTIF_HITS = args.output_file
TF_NAME_FILE = args.tf_name_file
DOMAIN_BED = args.domain_bed
RE_BED = args.RE_bed
CHIP_BED = args.CHIP_bed
GEOS_META = args.GEOS_meta

####################################------------CHECKPOINT 1-----------------###################################
# From here on the script will process the motif file that was produced from MOODS or provided as the -o option.
################################################################################################################
# The TF file must look something ilke the following with a header column that matches the following column IDs.
#
# Motif_ID		TF_Name
# M00116_1.97d	TFAP2B

print ("TACMAN: Parsing MOODS")

# Manipulate the list to get the selected feature
# create a datasframe from the TF_Name list
TF_df = pd.read_csv(TF_NAME_FILE, sep="\t", header=0)

# Crete a dictionary of motif names and tf names
dicted = snippets.dict_TF_df(TF_df)

a = snippets.parse_moods(MOTIF_HITS, dicted)

MOTIF_HITS_BED = MOTIF_HITS.replace(".txt", ".bed")
fish.animate()
bird.animate()
a.sort_values(by=['chr', "TF_start", "TF_end"], inplace=True)

a.to_csv(MOTIF_HITS_BED, sep="\t", index=False, header=False)

####################################------------CHECKPOINT 3-----------------###################################
# From here the parsed Moods file will be in a BED format and will be intersected with the RE lists of interest and also the domains.
################################################################################################################
print ("TACMAN: Finding Intersections")
print ("################################################################################################################")

#TF_domain_intersect = snippets.bed_intersect(MOTIF_HITS_BED, DOMAIN_BED, False)

#TF_domain_RE_intersect_df =  snippets.bed_intersect(RE_BED, TF_domain_intersect, True)

#TF_domain_RE_intersect_df["RE_ID"] = TF_domain_RE_intersect_df[0] + "_" + TF_domain_RE_intersect_df[1].apply(str) + "_" + TF_domain_RE_intersect_df[2].apply(str)
#TF_domain_RE_intersect_df["SUB_ID"] = TF_domain_RE_intersect_df[0] + "_" + TF_domain_RE_intersect_df[18].apply(str) + "_" + TF_domain_RE_intersect_df[19].apply(str)

#TF_domain_RE_intersect_df.drop([1, 2, 5, 12, 13, 14, 16, 17, 18, 19, 20], axis=1, inplace=True)

#TF_domain_RE_intersect_df = TF_domain_RE_intersect_df[[0, 6, 7, 8, 9, 10, 11, 15, "RE_ID", "SUB_ID", 3, 4]]

#TF_domain_RE_intersect_df.to_csv("TF_RE_SUB.bed", sep="\t", header = False, index=False)

MOTIF_PATH = os.path.abspath("TF_RE_SUB.bed")
####################################------------CHECKPOINT 4-----------------###################################
# The next step is to import the ChIP files from MARIO and arrange them by run mode, collect them into groups of common TF_start
# then take the nth percentile and return files either merged or unmerged for each tfself.
# This will also include parsing the GEOS metadatafile from MARIOS and using it to reference.
################################################################################################################
print ("################################################################################################################")
print ("TACMAN: Parsing GEOS META")

geos_df = snippets.parse_GEOS(GEOS_META)

print ("TACMAN: Parsing ChIP files")

CHIP_group_df = snippets.parse_CHIP(CHIP_BED, geos_df)

print ("TACMAN: Making ChIP folder")

snippets.make_set_dir("chip")

print ("TACMAN: Parsing the 75th percentile")

#snippets.percentile_parse(CHIP_group_df, 4)

print ("TACMAN: Making BED folder")

snippets.make_set_dir("merged")

print ("TACMAN: Merging duplicate BEDS")

#snippets.merge_replicate_BEDS(CHIP_group_df, MOTIF_PATH)

print ("TACMAN: Files merged by mode")

####################################------------CHECKPOINT 5-----------------###################################
# The next step is to take the merged BED files and then create a composite BED file to BIN based on MODE
################################################################################################################
print ("TACMAN: Binning Genome Based on MODES")
print ("TACMAN: Building reference files")

MODES = ["MODE1", "MODE2", "MODE3", "MODE4"]
MODES = pd.DataFrame(MODES)

files_list = []

#In this line of code you will break down the different MODES into sets that will be binned. The next step is sorting.
MODES[0].apply(snippets.bin_group_collect, args = (CHIP_group_df, MOTIF_PATH, files_list))

snippets.bin_groups_write_files(files_list)

####################################------------CHECKPOINT 6-----------------###################################
# The next step is to take the merged BED files and then create a composite BED file to BIN based on MODE
################################################################################################################
print ("TACMAN: Intersecting binned files with ChIP files and MOODS predictions")

snippets.intersect_bin(MOTIF_PATH)
