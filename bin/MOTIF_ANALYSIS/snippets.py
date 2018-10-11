#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import argparse
import pandas as pd
import seaborn as sns
import pybedtools
import math
import numpy as np
import scipy
import glob
import colored
from itertools import groupby, chain
from colored import stylize
from collections import OrderedDict

MODES = ["MODE1", "MODE2", "MODE3", "MODE4"]
MODES = pd.DataFrame(MODES)

#COLOR PARAMETERS
moods_color = colored.fg(35) + colored.attr(1)
BED_INTERSECT_color = colored.fg(141) + colored.attr(1)
DOMAIN_PARSE_color = colored.fg(202) + colored.attr(1)
MARIO_PARSE_color = colored.fg(14) + colored.attr(1)
MAKE_DIR_color = colored.fg(201) + colored.attr(1)
MERGE_BED_color = colored.fg(196) + colored.attr(1)
PARSE_PERCENTILE_color = colored.fg(15) + colored.attr(1)

def bin_parse(x):
    BIN_FILE = []
    peaks = x["ID"]
    bins = x['bins']
    start = x['Start']
    chrom = str(x['Chr'])
    bins_zero = bins - 1
    counter = 0

    print ("Binning BED: Writing Bins for: " + peaks)
    for i in range(int(bins_zero)):
        bin_track = counter * 100
        bin_start = start + bin_track
        bin_end = bin_start + 99

        tmp3 = chrom, bin_start, bin_end, peaks, bins

        BIN_FILE.append(tmp3)

        counter = counter + 1

    return BIN_FILE

def bed_intersect(x, y, z):
    print (stylize("BED Intersect: Importing BED1", BED_INTERSECT_color))
    a = pybedtools.BedTool(x)

    print (stylize("BED Intersect: Importing BED2", BED_INTERSECT_color))
    b = pybedtools.BedTool(y)

    print (stylize("BED Intersect: Intersecting BEDS", BED_INTERSECT_color))
    a_and_b = a.intersect(b, wa=True, wb=True)

    print (stylize("BED Intersect: Writing BEDS", BED_INTERSECT_color))
    c = a_and_b.moveto("TF_domain.bed")

    print (stylize("BED Intersect: Converting to frame", BED_INTERSECT_color))
    df = c.to_dataframe()

    if z == True:
        return df
    else:
        return c

def bin_group_collect(x, df, path, files_list):
    f_name = "DHS_CHIP_" + x + ".txt"
    files_list.append(f_name)
    frame = pd.DataFrame()
    list_df = []

    #create a dataframe of the TFs that have the same MODE
    print ("Bin Rerence Bed: Framing " + x + " BED files for processing")
    mode_df = df[df.MODE == x]

    # create a column with the name of the bed file for each MODE
    print ("Bin Rerence Bed: Extracting extensions")
    mode_df["bed_name"] = mode_df["TF_Name"] + "_" + mode_df["MODE"] + ".bed"

    # Create a list of bed files the concactenate
    print ("Bin Rerence Bed: Extracting unique BED files")
    bed_list = mode_df["bed_name"].unique().tolist()
    bed_list.append(path)

    for fname in bed_list:
        print ("Bin Rerence Bed: Creating DF of " + fname)
        bin_df = pd.read_csv(fname, sep="\t", usecols=[0, 1, 2], header=None)
        list_df.append(bin_df)

    print ("Bin Rerence Bed: Concatenating BED Dataframes")
    frame = pd.concat(list_df)

    print ("Bin Rerence Bed: Sorting: " + f_name)
    frame.sort_values(by=[0, 1, 2], inplace=True)

    print ("Bin Rerence Bed: Merging overlaps: " + f_name)
    a = pybedtools.BedTool.from_dataframe(frame)
    c = a.merge()
    df = c.moveto(f_name)

    print ("Bin Rerence Bed: Creating DF of: " + f_name)
    df = pd.read_csv(f_name, sep="\t", header=None)

    df.columns = ["Chr", "Start", "Stop"]

    print ("Bin Rerence Bed: Calculating BIN length")
    df["len"] = df["Stop"] - df["Start"] + 1

    print ("Bin Rerence Bed: Calculating BIN number")
    df["bins"] = df["len"]/100

    print ("Bin Rerence Bed: Applying ceilings to BIN numbers")
    df["bins"] = df["bins"].apply(math.ceil)
    df["ID"] = df["Chr"] + "_" + df["Start"].apply(str) + "_" + df["Stop"].apply(str)

    name = f_name + "_BINS.bed"

    print ("Binning BED: Working on: " + f_name)

    BIN_FILE = df.apply(a_bin_parse, axis=1)
    print ("Binning BED: Writing : " + f_name)

    BIN_FILE.to_csv(name, sep="\t", header=False, index=False)

def domain_parse(x, y):
        print (stylize("Domain Parsing: Extracting domain ID", DOMAIN_PARSE_color))
        x["SUB_ID"] = x[5] + "_" + x[6].apply(str) + "_" + x[7].apply(str)

        print (stylize("Domain Parsing: Ordering frame", DOMAIN_PARSE_color))
    	x = x[[0, 1, 2, 4, 3, 8, "SUB_ID"]]

        print (stylize("Domain Parsing: Replacing extensions", DOMAIN_PARSE_color))
    	MOODS_HITS_FN = y.replace(".txt", ".bed")

        print (stylize("Domain Parsing: Writing to CSV", DOMAIN_PARSE_color))
    	x.to_csv(MOODS_HITS_FN, sep="\t", header = False, index=False)

        return MOODS_HITS_FN

def intersect_bin(x):
    for m in MODES:
        f_name = "DHS_CHIP_" + m + ".txt_BINS.bed"

        y = glob.glob("*" + m + ".bed")

        print ("Intersect BED: List to intersect: " + y)

        a = pybedtools.BedTool(f_name)
        b = pybedtools.BedTool(y)

        print ("Intersect BED: Working on: " + f_name)

        a_and_b = a.intersect(b, wa=True, wb=True)

        c = a_and_b.moveto("CHIP_DHS_RE_SUB_BINS_" + m + ".bed")

        d = pybedtools.BedTool(x)

        print ("Intersect BED: Working on: CHIP")

        a_and_d = a.intersect(d, wa=True, wb=True)

        e = a_and_d.moveto("MOODS_DHS_RE_SUB_BINS_" + m + ".bed")

def make_set_dir(x):
    cwd = os.getcwd()

    dir = cwd + "/" + x
    print (stylize("Making Directory: " + dir, MAKE_DIR_color))

    if not os.path.exists(dir):
            os.makedirs(dir)

    print (stylize("Making Directory: Changing working directory to: " + dir, MAKE_DIR_color))
    os.chdir(dir)

def parse_MARIO(x, y):
    print (stylize("MARIO Parsing: Reading in GEOS Meta File", MARIO_PARSE_color))
    df = pd.read_csv(x, sep="\t", header=None, usecols=[1, 13], names=["Sample_Name", "info"])
    CHIP_df = pd.DataFrame.from_dict(y)

    print (stylize("MARIO Parsing: Splitting strings", MARIO_PARSE_color))
    df["TF"] = df["info"].str.split(":").str[3]
    df["Type"] = df["info"].str.split(":").str[1]

    print (stylize("MARIO Parsing: Dropping non-TF ChIP Files", MARIO_PARSE_color))

    histone_mod_list = [  "H3K9me3", "H3K79me2"]

    df = df[df.Type != "viral_protein"]
    df = df[df.TF != "H3K27ac"]
    df = df[df.TF != "H3K27me3"]
    df = df[df.TF != "H3K36me3"]
    df = df[df.TF != "H3K4me1"]
    df = df[df.TF != "H3K4me2"]
    df = df[df.TF != "H3K4me3"]
    df = df[df.TF != "H3K9ac"]
    df = df[df.TF != "H3K9me3"]
    df = df[df.TF != "H3K79me2"]

    print (stylize("MARIO Parsing: Creating DF ChIP name associations based on GEOS Meta File", MARIO_PARSE_color))
    CHIP_df["Basename"] = CHIP_df[0].apply(os.path.basename)
    CHIP_df["File_ext"] = CHIP_df.Basename.str.split("_").str[-1]
    CHIP_df["MODE"] = CHIP_df.File_ext.str.split(".").str[0]
    CHIP_df["Sample_Name"] = CHIP_df.Basename.str.split("_").str[0]

    print (stylize("MARIO Parsing: Droping unused columns", MARIO_PARSE_color))
    CHIP_df.drop("File_ext", axis=1, inplace=True)
    CHIP_df.drop("Basename", axis=1, inplace=True)

    print (stylize("MARIO Parsing: Mapping Names", MARIO_PARSE_color))
    dicted = dict(zip(df.Sample_Name, df.TF))

    CHIP_df["TF_Name"] = CHIP_df["Sample_Name"].map(dicted)
    CHIP_df.dropna(inplace=True)
    CHIP_df["file_path"] = CHIP_df[0]

    return CHIP_df

def parse_moods(x, TF_NAME_FILE, ofd):

    def dict_TF_df(x):
        df = pd.DataFrame()

        df["Motif"] = x["Motif_ID"]
        df["TF"] = x["TF_Name"]

        df.set_index("Motif")

        dicted = dict(zip(df.Motif, df.TF))

        return dicted

    # create a datasframe from the TF_Name list
    TF_df = pd.read_csv(TF_NAME_FILE, sep="\t", header=0)

    # Crete a dictionary of motif names and tf names
    dicted = dict_TF_df(TF_df)

    column_names = ['PEAK_ID', 'PWM_FILE', 'TF_START', 'STRAND', 'MATCH_SCORE', 'MOTIF_SEQ']
    column_types = {
        'PEAK_ID': "category",
        'PWM_FILE': "category",
        'TF_START': "uint16",
        'STRAND': "category",
        'MATCH_SCORE': "float32",
        'MOTIF_SEQ': "category"
    }

    print (stylize("Parsing MOODS: Reading in CSV", moods_color))
    df = pd.read_csv(x, usecols=[0,1,2,3,4,5], names=column_names, dtype=column_types, header=None, sep="|")

    print (stylize("Parsing MOODS: Replacing some extensions", moods_color))
    df['PWM_FILE'] = df['PWM_FILE'].str.replace('_JASPAR.txt.pfm', '')

    print (stylize("Parsing MOODS: Splitting and creating new frame", moods_color))
    df_tmp1 = df['PEAK_ID'].str.split(":", expand=True)

    print (stylize("Parsing MOODS: Splitting and creating new frame2", moods_color))
    df_tmp2 = df_tmp1[1].str.split("-", expand=True)

    print (stylize("Parsing MOODS: Dropping another column", moods_color))
    df.drop(columns=['PEAK_ID'], inplace=True)

    print (stylize("Parsing MOODS: Setting Values", moods_color))
    df["chr"] = df_tmp1[0]
    df["start"] = df_tmp2[0]
    df["stop"] = df_tmp2[1]

    print (stylize("Parsing MOODS: Deleting unused frames", moods_color))
    del df_tmp1
    del df_tmp2

    print (stylize("Parsing MOODS: Calculating TF motif start positions", moods_color))
    df["TF_start"] = df["start"].apply(int) + 1 + df["TF_START"]

    print (stylize("Parsing MOODS: Calculating TF motif end positions", moods_color))
    df["TF_end"] = df["TF_start"] + df["MOTIF_SEQ"].str.len() - 1

    print (stylize("Parsing MOODS: Extracting Peak ID", moods_color))
    df["PEAK_ID"] = df["chr"] + "_" + df["start"].apply(str) + "_" + df["stop"].apply(str)

    print (stylize("Parsing MOODS: Extracting Motif ID", moods_color))
    df["MOTIF_POS"] = df["chr"] + "_" + df["TF_start"].apply(str) + "_" + df["TF_end"].apply(str)

    print (stylize("Parsing MOODS: Calculating motif length", moods_color))
    df["MOTIF_LEN"] = df["TF_end"] - df["TF_start"] + 1

    print (stylize("Parsing MOODS: Mapping TF names", moods_color))
    df["TF_Name"] = df["PWM_FILE"].map(dicted)

    print (stylize("Parsing MOODS: Dropping unused columns", moods_color))
    #df.drop(columns=["TF_START", "PWM_FILE", "MATCH_SCORE", "MOTIF_POS", "MOTIF_LEN", "MOTIF_SEQ", "STRAND", "start", "stop"], inplace=True)

    df = df[["chr", "TF_start", "TF_end", "PEAK_ID", "TF_Name"]]

    MOODS_HITS_BED = os.path.basename(x)
    MOODS_HITS_BED = MOODS_HITS_BED.replace(".txt", ".bed")

    print (stylize("Parsing MOODS: Sorting columns", moods_color))
    df.sort_values(by=['chr', "TF_start", "TF_end"], inplace=True)

    print (stylize("Parsing MOODS: Writing CSV", moods_color))
    df.to_csv(MOODS_HITS_BED, sep="\t", index=False, header=False)

    return MOODS_HITS_BED

def parse_percentile(x, y):
    base = os.path.basename(x)
    print (stylize("Parse 75th percentile from: " + str(base), PARSE_PERCENTILE_color))

    df = pd.read_csv(str(x), sep="\t", header=None)
    col_num = len(df.columns)

    if col_num == 11:
        per_len = (len(df[7]))/4

        df.sort_values(by=[7], ascending=False, inplace=True)
        df = df.head(per_len)

        df.to_csv(str(base))

    elif col_num == 13:
        per_len = (len(df[12]))/4

        df.sort_values(by=[12], ascending=False, inplace=True)
        df = df.head(per_len)

        df.to_csv(str(base))

def parse_replicate_BEDS(x, y):
    def parse_percentile(p):
        print (stylize("Parse 75th percentile from: " + str(base), PARSE_PERCENTILE_color))
        df = pd.read_csv(str(p), sep="\t", header=None)
        col_num = len(df.columns)

        if col_num == 11:
            per_len = (len(df[7]))/4
            df.sort_values(by=[7], ascending=False, inplace=True)
            df = df.head(per_len)

            return df

        elif col_num == 13:
            per_len = (len(df[12]))/4
            df.sort_values(by=[12], ascending=False, inplace=True)
            df = df.head(per_len)

            return df

    def parse_replicates(i):
        tf_df = x[x.TF_Name == i]
        print (stylize("Merge BEDS: Merging " + i, MERGE_BED_color))

        for j in MODES[0]:
            mode_df = tf_df[tf_df.MODE == j]

            tf_filename = i + "_" + str(j) + ".bed"
            with open(tf_filename, 'w') as outfile:
                for fname in mode_df["file_path"]:
                    with open(fname) as infile:
                        outfile.write(infile.read())

            df_sort = pd.read_csv(tf_filename, sep="\t", header=None, usecols=[0,1,2,3])
            df_sort.sort_values(by=[0, 1, 2], inplace=True)
            df_sort.to_csv(tf_filename, sep="\t", index=False, header=False)

            a = pybedtools.BedTool(tf_filename)
            c = a.merge()
            d = c.moveto(tf_filename)

    def parse_single(k):
        tf_df = x[x.TF_Name == k]
        print (stylize("Merge BEDS: Parsing TF single file " + k, MERGE_BED_color))

        for l in MODES[0]:
            mode_df = tf_df[tf_df.MODE == l]

            tf_filename = k + "_" + l + ".bed"
            with open(tf_filename, 'w') as outfile:
                for fname in mode_df["file_path"]:
                    with open(fname) as infile:
                        outfile.write(infile.read())

    base = os.path.basename(p)

    print (stylize("Merge BEDS: Creating reference frames", MERGE_BED_color))
    TF_file_counts = x.groupby(["TF_Name", "MODE"]).count()
    TF_file_counts.drop([0], axis=1, inplace=True)
    TF_file_counts.reset_index(inplace=True)

    TF_replicates = TF_file_counts[TF_file_counts.Sample_Name > 1].unique()
    TF_singles = TF_file_counts[TF_file_counts.Sample_Name == 1].unique()

    print (TF_replicates)
    print (TF_singles)

    TF_replicates["TF_Name"].apply(parse_replicates)
    TF_singles["TF_Name"].apply(parse_single)
