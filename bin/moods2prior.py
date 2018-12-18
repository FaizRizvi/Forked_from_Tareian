#!/usr/bin/env python
from __future__ import print_function

import pandas as pd
import seaborn as sns
import os
import sys
import argparse
import pybedtools
import numpy as np

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

# Create the dataframe that will be used to get the dictionary of TF names
df_TF = pd.read_csv(args.IN_TF, sep="\t", header=0)

# Get the base names and set as the index ---- This might need to be changed later
df_TF_names = pd.DataFrame()
df_TF_names["Motif"] = df_TF["Motif_ID"].str.replace("_1.97d", "")
df_TF_names["TF"] = df_TF["TF_Name"]
df_TF_names.set_index("Motif")

dict_TF = dict(zip(df_TF_names.Motif, df_TF_names.TF))

del df_TF
del df_TF_names

##########################-----------LOAD GENE INFORMATION------------##############################
print("Loading and Parsing Gene Information")

# Create the dataframe that will be used to get the dictionary of Gene names
df_gene = pd.read_csv(args.GENE, sep="\t", header=0)

# Get the base names and set as the index ---- This might need to be changed later
df_gene_names = pd.DataFrame()
df_gene_names["UCSC"] = df_gene["#mm10.knownToEnsembl.name"]
df_gene_names["Symbol"] = df_gene["mm10.kgXref.geneSymbol"]
df_gene_names.set_index("UCSC")

dict_gene = dict(zip(df_gene_names.UCSC, df_gene_names.Symbol))

del df_gene
del df_gene_names

##########################-----------LOAD Domain Information------------##############################
print("Loading and Parsing Domain Information")

bed_TSS = pybedtools.BedTool(args.TSS)
bed_domain = pybedtools.BedTool(args.DOMAIN)
bed_TSS_and_bed_domain = bed_TSS.window(bed_domain, w=1)

del bed_domain
del bed_TSS

#write_bed = bed_TSS_and_bed_domain.moveto(name_domain)

##########################-----------LOAD MOODS DATA------------##############################
print("Loading and Parsing MOODS Information")

column_names = ['ID', 'Motif', 'TF_pos', 'Strand', 'Match_score', 'Motif_sequence']

column_types = {
    'ID': "category",
    'Motif': "category",
    'TF_pos': "uint16",
    'Strand': "category",
    'Match_score': "float32",
    'Motif_sequence': "category"
}

# Create a DF that will read in the moods input
df_moods = pd.read_csv(args.IN_MOODS, sep="|", usecols=[0,1,2,3,4,5], header=None, names=column_names, dtype=column_types)

# Create a temp df that will house the expanded MOTIF information
tmp1 = df_moods["Motif"].str.split('_', 1, expand=True)

df_moods.drop("Motif", axis=1)
df_moods["Motif"] = tmp1[0]

del tmp1

df_moods["TF_Name"] = df_moods["Motif"].map(dict_TF)

# Create a temp df that will house the expanded SEQUENCE information
tmp2 = df_moods["ID"].str.split(':', 1, expand=True)
tmp3 = tmp2[1].str.split("-", 1, expand=True)

df_moods["Chr"] = tmp2[0]
df_moods["Start"] = tmp3[0]
df_moods["Stop"] = tmp3[1]

del tmp2
del tmp3

df_moods["TF_start"] = df_moods["Start"].apply(int) + 1 + df_moods["TF_pos"].apply(int)
df_moods["TF_end"] = df_moods["TF_start"] + df_moods["Motif_sequence"].str.len() - 1
df_moods["PEAK_ID"] = df_moods["Chr"] + "_" + df_moods["Start"].map(str) + "_" + df_moods["Stop"].map(str)

order= ["Chr", "TF_start", "TF_end", "TF_Name", "Strand", "PEAK_ID"]

df_moods = df_moods[order]

bed_moods = pybedtools.BedTool.from_dataframe(df_moods)
bed_TSS_and_moods = bed_TSS_and_bed_domain.window(bed_moods, w=10000)

del df_moods
del bed_moods

df_intersect = bed_TSS_and_moods.to_dataframe(names=["TSS_chr", "TSS_start", "TSS_end", "Gene_name", "Strand", "Domain_chr", "Domain_start", "Domain_end", "Domain", "Nope", "dot", "s_c", "s_s", "s_e", "RGB", "TF_chr", "TF_start", "TF_end", "TF_Motif", "TF_strand", "Peak_ID"])

#d = bed_TSS_and_moods.moveto(name_TSS)
del bed_TSS_and_moods

df_intersect["SUB_ID"] = df_intersect["Domain_chr"] + "_" + df_intersect["Domain_start"].apply(str) + "_" + df_intersect["Domain_end"].apply(str) + "_" + df_intersect["Domain"]

keep = ["Gene_name", "TF_end", "SUB_ID", "TSS_chr"]

df_intersect = df_intersect[keep]

df_intersect["Gene_Symbol"] = df_intersect["Gene_name"].map(dict_gene)

##########################-----------Count and Pivot by distance------------##############################
print("Counting and pivoting DF")

df_gb = df_intersect.groupby(["Gene_Symbol", "TF_end"])["TSS_chr"].count()

df_gb.reset_index(inplace=True)

df_piv = df_gb.pivot_table(index="Gene_Symbol", columns="TF_end", values="TSS_chr", aggfunc=np.sum)

df_piv = df_piv.fillna(0)

df_piv_binary = df_piv.apply(lambda x: [y if y <= 1 else 1 for y in x])

df_piv.to_csv(name_piv, sep="\t")
df_piv_binary.to_csv(name_binary, sep="\t")

del df_gb
del df_piv

##########################-----------Count and Pivot by TAD------------##############################
print("Counting and pivoting DF")

df_gb_tad = df_intersect.groupby(["SUB_ID", "Gene_Symbol", "TF_end"])["TSS_chr"].count()

df_gb_tad.reset_index(inplace=True)

df_piv_tad = df_gb_tad.pivot_table(index="Gene_Symbol", columns="TF_end", values="TSS_chr", aggfunc=np.sum)

df_piv_tad = df_piv.fillna(0)

df_piv_binary_tad = df_piv_tad.apply(lambda x: [y if y <= 1 else 1 for y in x])

df_piv_tad.to_csv(name_piv, sep="\t")
df_piv_binary_tad.to_csv(name_binary, sep="\t")