#!/usr/bin/env python
from __future__ import print_function

import pandas as pd
import seaborn as sns
import os
import sys
import argparse
import pybedtools
import numpy as np
import snippets

##########################-----------ARGUMENT PARSER------------##############################
"""Set up the argument parser with all of the options and defaults set up."""
parser = argparse.ArgumentParser()

parser.add_argument("-m", "--MOODS", dest='IN_MOODS', help="Input MOODS file: sep = |", required=True)
parser.add_argument("-t", "--TF", dest='IN_TF', help="Input TF information file: sep = \t", required=True)
parser.add_argument("-T", "--TSS", dest='TSS', help="Input TSS information file: sep = \t", required=True)
parser.add_argument("-g", "--GENE", dest='GENE', help="Input GENE information file: sep = \t", required=True)
parser.add_argument("-o", "--output", dest='OUT_NAME', help="Output file name", required=True)
parser.add_argument("-d", "--DOMAIN", dest='DOMAIN', help="Domain List", required=True)

# set the arguments from the command line to variables in the args object
args = parser.parse_args()

name_binary = args.OUT_NAME + "_binary.txt"
name_TSS = args.OUT_NAME + "_TSS_10kb.bed"
name_domain = args.OUT_NAME + "_domain.txt"
name_piv = args.OUT_NAME + "_prior_counts.txt"

##########################-----------LOAD TF INFORMATION------------##############################
print("Loading and Parsing TF Information")

dict_TF = TFmeta2dict(args.IN_TF)

##########################-----------LOAD GENE INFORMATION------------##############################
print("Loading and Parsing Gene Information")

dict_gene = genemeta2dict(args.GENE)

##########################-----------LOAD Domain Information------------##############################
print("Loading and Parsing Domain Information")

bed_TSS_and_bed_domain = pybedtools_window(args.TSS, args.DOMAIN, 1, False, name_domain)

##########################-----------LOAD MOODS DATA------------##############################
print("Loading and Parsing MOODS Information")

df_moods = parse_moods(args.IN_MOODS)

bed_moods = pybedtools.BedTool.from_dataframe(df_moods)
bed_TSS_and_moods = bed_TSS_and_bed_domain.window(bed_moods, w=10000)

df_intersect = bed_TSS_and_moods.to_dataframe(names=["TSS_chr", "TSS_start", "TSS_end", "Gene_name", "Strand", "Domain_chr", "Domain_start", "Domain_end", "Domain", "Nope", "dot", "s_c", "s_s", "s_e", "RGB", "TF_chr", "TF_start", "TF_end", "TF_Motif", "TF_strand", "Peak_ID"])

del df_moods
del bed_moods
del bed_TSS_and_moods

df_intersect["SUB_ID"] = df_intersect["Domain_chr"] + "_" + df_intersect["Domain_start"].apply(str) + "_" + df_intersect["Domain_end"].apply(str) + "_" + df_intersect["Domain"]

keep = ["Gene_name", "TF_end", "SUB_ID", "TSS_chr"]

df_intersect = df_intersect[keep]

df_intersect["Gene_Symbol"] = df_intersect["Gene_name"].map(dict_gene)

##########################-----------Count and Pivot by distance------------##############################
print("Counting and pivoting DF")

snippets.groupby_pivot(df_intersect, ["Gene_Symbol", "TF_end"], "TSS_chr")

##########################-----------Count and Pivot by TAD------------##############################
print("Counting and pivoting DF")

snippets.groupby_pivot(df_intersect, ["SUB_ID", "Gene_Symbol", "TF_end"], "TSS_chr")
