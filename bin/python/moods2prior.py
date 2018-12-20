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
parser.add_argument("-d", "--DOMAIN", dest='DOMAIN', help="Input TAD information file: sep = \t", required=True)

# set the arguments from the command line to variables in the args object
args = parser.parse_args()

name_TSS = args.OUT_NAME + "_TSS_10kb.bed"
name_domain = args.OUT_NAME + "_domain.txt"
name_piv = args.OUT_NAME + "_prior_sum.txt"
name_binary = args.OUT_NAME + "_prior_binary.txt"
name_TAD = args.OUT_NAME + "_prior_TAD_sum.txt"
name_TAD_binary = args.OUT_NAME + "_prior_TAD_binary.txt"

##########################-----------LOAD TF INFORMATION------------##############################
print("Loading and Parsing TF Information")

dict_TF = snippets.TFmeta2dict(args.IN_TF)

##########################-----------LOAD GENE INFORMATION------------##############################
print("Loading and Parsing Gene Information")

dict_gene = snippets.genemeta2dict(args.GENE)

##########################-----------LOAD MOODS DATA------------##############################
print("Loading and Parsing MOODS Information")

bed_moods = snippets.parse_moods(args.IN_MOODS, dict_TF)

##########################-----------Intersect 10kb of TSS------------##############################
print("Intersecting BEDS")

bed_TSS_and_bed_moods = snippets.pybedtools_window(args.TSS, bed_moods, 10000, True)
df_intersect = bed_TSS_and_bed_moods.to_dataframe(names=["Gene_chr", "Gene_start", "Gene_end", "Gene_name", "Gene_strand", "TF_chr", "TF_start", "TF_end", "TF_name", "PEAK_ID"])

keep = ["Gene_name", "TF_name", "PEAK_ID"]
df_intersect = df_intersect[keep]

df_intersect["Gene_Symbol"] = df_intersect["Gene_name"].map(dict_gene)

##########################-----------Intersect TADS------------##############################
print("Intersecting BEDS")

bed_domain = snippets.pybedtools_intersect(args.TSS, args.DOMAIN)

bed_domain_moods = bed_moods.intersect(bed_domain, wa=True, wb=True)

df_domain_mood = bed_domain_moods.to_dataframe(names=["TF_chr", "TF_start", "TF_end", "TF_name", "PEAK_ID", "Gene_chr", "Gene_start", "Gene_end", "Gene_name", "Gene_strand", "Dom_chr", "Dom_start", "Dom_end", "Dom_Name", "Num", "Dot", "dom_start", "dom_end", "rgb"])

df_domain_mood["DOMAIN_ID"] = df_domain_mood["Dom_chr"] + "_" + df_domain_mood["Dom_start"].apply(str) + "_" + df_domain_mood["Dom_end"].apply(str) + "_" + df_domain_mood["Dom_Name"]

df_domain_mood["Gene_Symbol"] = df_domain_mood["Gene_name"].map(dict_gene)

keep2 = ["TF_name", "Gene_Symbol", "PEAK_ID", "DOMAIN_ID"]

df_domain_mood = df_domain_mood[keep2]

##########################-----------Count and Pivot------------##############################
print("Counting and pivoting DF")

snippets.groupby_pivot(df_intersect, ["Gene_Symbol", "TF_name"], name_piv, name_binary)

snippets.groupby_pivot(df_domain_mood, ["DOMAIN_ID", "Gene_Symbol", "TF_name"], name_TAD, name_TAD_binary)
