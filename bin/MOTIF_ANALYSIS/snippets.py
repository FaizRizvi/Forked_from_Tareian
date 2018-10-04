import MOODS.scan
import MOODS.tools
import MOODS.parsers
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
import time
from tqdm import tqdm

from itertools import groupby, chain

MODES = ["MODE1", "MODE2", "MODE3", "MODE4"]
MODES = pd.DataFrame(MODES)

def a_bin_parse(x):
    BIN_FILE = []
    peaks = x["ID"]
    bins = x['bins']
    start = x['Start']
    chrom = str(x['Chr'])
    bins_zero = bins - 1
    counter = 0

    print ("Binning BED: Writing Bins for: " + peaks)
    for i in tqdm(range(int(bins_zero))):
        bin_track = counter * 100
        bin_start = start + bin_track
        bin_end = bin_start + 99

        tmp3 = chrom, bin_start, bin_end, peaks, bins

        BIN_FILE.append(tmp3)

        counter = counter + 1

    return BIN_FILE

def bed_intersect(x, y, z):
    print ("BED Intersect: Importing BED1")
    a = pybedtools.BedTool(x)

    print ("BED Intersect: Importing BED2")
    b = pybedtools.BedTool(y)

    print ("BED Intersect: Intersecting BEDS")
    a_and_b = a.intersect(b, wa=True, wb=True)

    print ("BED Intersect: Writing BEDS")
    c = a_and_b.moveto("TF_domain.bed")

    print ("BED Intersect: Creating BED files")
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

def dict_TF_df(x):
    df = pd.DataFrame()

    df["Motif"] = x["Motif_ID"]
    df["TF"] = x["TF_Name"]

    df.set_index("Motif")

    dicted = dict(zip(df.Motif, df.TF))

    return dicted

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

    if not os.path.exists(dir):
            os.makedirs(dir)

    os.chdir(dir)

def merge_replicate_BEDS(x, y):
    replicate_groups = x["TF_Name"].unique()

    rep_group_counts = x.groupby(["TF_Name", "MODE"]).count()
    rep_group_counts.drop([0], axis=1, inplace=True)
    rep_group_counts.reset_index(inplace=True)

    rep_group_replicates = rep_group_counts[rep_group_counts.Sample_Name > 1]
    rep_group_singles = rep_group_counts[rep_group_counts.Sample_Name == 1]
    rep_group_all = rep_group_counts[rep_group_counts.Sample_Name >= 1]

    unique_TF_Rep = rep_group_replicates["TF_Name"].unique()
    unique_TF_single = rep_group_singles["TF_Name"].unique()
    unique_TF_all = rep_group_all["TF_Name"].unique()

    for i in unique_TF_Rep:
        tf_df = x[x.TF_Name == i]

        for j in MODES:
            mode_df = tf_df[tf_df.MODE == j]

            mode_paths = mode_df["file_path"]

            tf_filename = i + "_" + j + ".bed"

            print ("Merging " + tf_filename)

            with open(tf_filename, 'w') as outfile:
                for fname in mode_paths:
                    with open(fname) as infile:
                        outfile.write(infile.read())

            df_sort = pd.read_csv(tf_filename, sep="\t", header=None, usecols=[0,1,2,3])
            df_sort = df_sort[[0,1,2,3]]
            df_sort.sort_values(by=[0, 1, 2], inplace=True)

            df_sort.to_csv(tf_filename, sep="\t", index=False, header=False)

            a = pybedtools.BedTool(tf_filename)
            c = a.merge()
            d = c.moveto(tf_filename)

    for k in unique_TF_single:
        tf_df = x[x.TF_Name == k]

        for l in MODES:
            mode_df = tf_df[tf_df.MODE == l]

            mode_paths = mode_df["file_path"]

            tf_filename = k + "_" + l + ".bed"
            print ("Parsing TF single file " + tf_filename)
            with open(tf_filename, 'w') as outfile:
                for fname in mode_paths:
                    with open(fname) as infile:
                        outfile.write(infile.read())

def parse_CHIP(x, y):
    df = pd.DataFrame.from_dict(x)

    df["Basename"] = df[0].apply(os.path.basename)
    df["File_ext"] = df.Basename.str.split("_").str[-1]
    df["MODE"] = df.File_ext.str.split(".").str[0]
    df["Sample_Name"] = df.Basename.str.split("_").str[0]

    df.drop("File_ext", axis=1, inplace=True)
    df.drop("Basename", axis=1, inplace=True)

    dicted = dict(zip(y.Sample_Name, y.TF))

    df["TF_Name"] = df["Sample_Name"].map(dicted)
    df.dropna(inplace=True)
    df["file_path"] = df[0]

    return df

def parse_GEOS(x):
    df = pd.read_csv(x, sep="\t", header=None, usecols=[1, 13], names=["Sample_Name", "info"])

    df["TF"] = df["info"].str.split(":").str[3]
    df["Type"] = df["info"].str.split(":").str[1]

    df.drop(["info"], axis=1, inplace=True)

    df = df[df.Type != "viral_protein"]
    histone_markers_list = ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3", "H3K79me2"]

    for i in histone_markers_list:
        df = df[df.TF != i]

    return df

def parse_moods(x, y):
    print ("Parsing MOODS: Reading in CSV")
    df = pd.read_csv(x, header=None, sep="|")

    print ("Parsing MOODS: Replacing some extensions")
    df[1] = df[1].str.replace('_JASPAR.txt.pfm', '')

    print ("Parsing MOODS: Dropping a column")
    df.drop(6, axis=1, inplace=True)

    print ("Parsing MOODS: Splitting and creating new frame")
    df_tmp1 = df[0].str.split(":", expand=True)

    print ("Parsing MOODS: Splitting and creating new frame2")
    df_tmp2 = df_tmp1[1].str.split("-", expand=True)

    print ("Parsing MOODS: Dropping another column")
    df.drop(columns=0, inplace=True)

    print ("Parsing MOODS: Setting Values")
    df["chr"] = df_tmp1[0]
    df["start"] = df_tmp2[0]
    df["stop"] = df_tmp2[1]

    print ("Parsing MOODS: Deleting unused frames")
    del df_tmp1
    del df_tmp2

    print ("Parsing MOODS: Setting column names")
    df.columns = ["MOTIF_ID", "TF_POS", "STRAND", "MATCH_SCORE", "MOTIF_SEQ", "chr", "start", "stop"]

    print ("Parsing MOODS: Calculating TF motif start positions")
    df["TF_start"] = df["start"].apply(int) + 1 + df["TF_POS"]

    print ("Parsing MOODS: Calculating TF motif end positions")
    df["TF_end"] = df["TF_start"] + df["MOTIF_SEQ"].str.len() - 1

    print ("Parsing MOODS: Extracting Peak ID")
    df["PEAK_ID"] = df["chr"] + "_" + df["start"].map(str) + "_" + df["stop"].map(str)

    print ("Parsing MOODS: Extracting Motif ID")
    df["MOTIF_POS"] = df["chr"] + "_" + df["TF_start"].map(str) + "_" + df["TF_end"].map(str)

    print ("Parsing MOODS: Calculating motif length")
    df["MOTIF_LEN"] = df["TF_end"] - df["TF_start"] + 1

    print ("Parsing MOODS: Mapping TF names")
    df["TF_Name"] = df["MOTIF_ID"].map(y)

    print ("Parsing MOODS: Dropping unused columns")
    df.drop(columns=["TF_POS", "MOTIF_ID", "MATCH_SCORE", "MOTIF_SEQ", "STRAND", "start", "stop"], inplace=True)

    return df

def percentile_parse(x, y):
    for i in x["file_path"]:
        df = pd.read_csv(str(i), sep="\t", header=None)
        col_num = len(df.columns)
        base = os.path.basename(i)

        if col_num == 11:
            per_len = (len(df[7]))/y

            top_sorted = df[7].sort_values(ascending=False).head(per_len)

            df.sort_values(by=[7], ascending=False, inplace=True)

            df = df.head(per_len)

            df.to_csv(str(base))

        elif col_num == 13:
            per_len = (len(df[12]))/4

            top_sorted = df[12].sort_values(ascending=False).head(per_len)

            df.sort_values(by=[12], ascending=False, inplace=True)

            df = df.head(per_len)

            df.to_csv(str(base))
