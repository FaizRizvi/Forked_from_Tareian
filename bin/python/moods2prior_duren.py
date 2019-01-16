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
# Set up the argument parser with all of the options and defaults set up.
parser = argparse.ArgumentParser()

parser.add_argument("-m", "--MOODS", dest='IN_MOODS', help="Input MOODS file: sep = |", required=True)
parser.add_argument("-t", "--TF", dest='IN_TF', help="Input TF information file: sep = \t", required=True)
parser.add_argument("-T", "--TSS", dest='TSS', help="Input TSS information file: sep = \t", required=True)
parser.add_argument("-g", "--GENE", dest='GENE', help="Input GENE information file: sep = \t", required=True)
parser.add_argument("-d", "--DOMAIN", dest='DOMAIN', help="Input TAD information file: sep = \t", required=True)
parser.add_argument("-o", "--OUT_DIR", dest='OUT_DIR', help="Output dir", required=True)

# set the arguments from the command line to variables in the args object
args = parser.parse_args()

# if the output directory does not exists. Create it
if not os.path.exists(args.OUT_DIR):
        os.makedirs(args.OUT_DIR)

# Change output directory 
os.chdir(args.OUT_DIR)

# Get the basename of the files
base = os.path.basename(args.IN_MOODS)
basename = os.path.splitext(base)[0]

# Set up the file names based on the basename
name_piv = basename + "_TSS_10kb_prior_sum.txt"
name_binary = basename + "_TSS_10kb_prior_binary.txt"
name_TAD = basename + "_prior_TAD_sum.txt"
name_TAD_binary = basename + "_prior_TAD_binary.txt"

##########################-----------LOAD TF INFORMATION------------##############################
print("Loading and Parsing TF Information")

dict_TF = snippets.TFmeta2dict(args.IN_TF)

##########################-----------LOAD GENE INFORMATION------------##############################
print("Loading and Parsing Gene Information")

dict_gene = snippets.genemeta2dict(args.GENE)

##########################-----------LOAD MOODS DATA------------##############################
print("Loading and Parsing MOODS Information")

bed_moods = snippets.parse_moods(args.IN_MOODS, dict_TF)

##########################-----------LOAD TSS DATA------------##############################
print("Loading TSS file and converting to BED")

bed_TSS = pybedtools.BedTool(args.TSS)

##########################-----------Intersect 10kb of TSS------------##############################
print("Intersecting BEDS")

bed_TSS_and_bed_moods = bed_TSS.window(bed_moods, w=10000)

columns = ["Gene_chr", "Gene_start", "Gene_end", "Gene_name", "Gene_strand", "TF_chr", "TF_start", "TF_end", "TF_name", "PEAK_ID"]
keep = ["Gene_name", "TF_name", "PEAK_ID"]

df_intersect = bed_TSS_and_bed_moods.to_dataframe(names=columns, usecols=keep)

df_intersect["Gene_Symbol"] = df_intersect["Gene_name"].map(dict_gene)

##########################-----------Intersect TADS------------##############################
#print("Intersecting BEDS")

#bed_domain = snippets.pybedtools_intersect(args.TSS, args.DOMAIN)

##bed_domain_moods = bed_moods.intersect(bed_domain, wa=True, wb=True)

#columns2 = ["TF_chr", "TF_start", "TF_end", "TF_name", "PEAK_ID", "Gene_chr", "Gene_start", "Gene_end", "Gene_name", "Gene_strand", "Dom_chr", "Dom_start", "Dom_end", "Dom_Name", "Num", "Dot", "dom_start", "dom_end", "rgb"]
#keep2 = ["TF_name", "PEAK_ID", "Dom_chr", "Dom_start", "Dom_end", "Dom_Name", "Gene_name"]

#df_domain_mood = bed_domain_moods.to_dataframe(names=columns2, usecols=keep2)

#df_domain_mood["DOMAIN_ID"] = df_domain_mood["Dom_chr"] + "_" + df_domain_mood["Dom_start"].apply(str) + "_" + df_domain_mood["Dom_end"].apply(str) + "_" + df_domain_mood["Dom_Name"]

#df_domain_mood["Gene_Symbol"] = df_domain_mood["Gene_name"].map(dict_gene)

##########################-----------Count and Pivot------------##############################
print("Counting and pivoting DF")

# Create the 10kb prior
snippets.generate_prior(df_intersect, ["Gene_Symbol", "TF_name"], name_piv, args.OUT_DIR, name_binary)

# Create the TAD based prior
#snippets.generate_prior(df_domain_mood, ["DOMAIN_ID", "Gene_Symbol", "TF_name"], name_TAD, args.OUT_DIR, False)
#snippets.generate_prior(df_domain_mood, ["DOMAIN_ID", "Gene_Symbol", "TF_name"], name_TAD_binary,  args.OUT_DIR, True)
