import MOODS.scan
import MOODS.tools
import MOODS.parsers

import os
import sys
import argparse
from itertools import groupby, chain
import pandas as pd
import seaborn as sns
import pybedtools
import math
import numpy as np
import scipy
import glob

MODES = ["MODE1", "MODE2", "MODE3", "MODE4"]

def bed_intersect(x, y, z):
    a = pybedtools.BedTool(x)
    b = pybedtools.BedTool(y)

    a_and_b = a.intersect(b, wa=True, wb=True)

    c = a_and_b.moveto("TF_domain.bed")

    df = c.to_dataframe()

    if z == True:
        return df
    else:
        return c

def bin_BED(x):
    for i in x:
        df = pd.read_csv(i, sep="\t", header=None)
        df.columns = ["Chr", "Start", "Stop"]

        df["len"] = df["Stop"] - df["Start"] + 1
        df["bins"] = df["len"]/100
        df["bins"] = df["bins"].apply(math.ceil)
        df["ID"] = df["Chr"] + "_" + df["Start"].apply(str) + "_" + df["Stop"].apply(str)

        name = i + "_BINS.bed"

        with open(name, "w") as fout:
          for peaks in df["ID"]:
            bins = df.loc[df['ID'] == peaks, 'bins'].values[0]
            start = df.loc[df['ID'] == peaks, 'Start'].values[0]
            chrom = df.loc[df['ID'] == peaks, 'Chr'].values[0]
            bins_zero = int(bins - 1)
            counter = 0
            while counter <= bins_zero:
                bin_track = counter * 100
                bin_start = start + bin_track
                bin_end = bin_start + 99
                counter = counter + 1
                fout.write(chrom + "\t" + str(bin_start) + "\t" + str(bin_end) + "\t" + peaks +  "\t" + str(bins) + "\n")

def build_reference_BEDS(x, y):
    files_to_sort = []

    for m in MODES:
        mode_df = x[x.MODE == m]
        mode_df["bed_name"] = mode_df["TF_Name"] + "_" + mode_df["MODE"] + ".bed"

        bed_list = mode_df["bed_name"].unique().tolist()
        bed_list.append(y)

        f_name = "DHS_CHIP_" + m + ".txt"

        files_to_sort.append(f_name)

        frame = pd.DataFrame()
        list_ = []

        for fname in bed_list:
            bin_df = pd.read_csv(fname, sep="\t", usecols=[0, 1, 2], header=None)
            list_.append(bin_df)

        frame = pd.concat(list_)
        frame.to_csv(f_name, index=False, header=False)

    return files_to_sort

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

        a = pybedtools.BedTool(f_name)
        b = pybedtools.BedTool(y)

        a_and_b = a.intersect(b, wa=True, wb=True)

        c = a_and_b.moveto("CHIP_DHS_RE_SUB_BINS_" + m + ".bed")

        d = pybedtools.BedTool(x)
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
    df = pd.read_csv(x, header=None, sep="|")

    df[1] = df[1].str.replace('_JASPAR.txt.pfm', '')

    df.drop(6, axis=1, inplace=True)

    df_tmp1 = df[0].str.split(":", expand=True)
    df_tmp2 = df_tmp1[1].str.split("-", expand=True)

    df.drop(columns=0, inplace=True)

    df["chr"] = df_tmp1[0]
    df["start"] = df_tmp2[0]
    df["stop"] = df_tmp2[1]

    df.columns = ["MOTIF_ID", "TF_POS", "STRAND", "MATCH_SCORE", "MOTIF_SEQ", "chr", "start", "stop"]

    df["TF_start"] = df["start"].apply(int) + 1 + df["TF_POS"]
    df["TF_end"] = df["TF_start"] + df["MOTIF_SEQ"].str.len() - 1
    df["PEAK_ID"] = df["chr"] + "_" + df["start"].map(str) + "_" + df["stop"].map(str)
    df["MOTIF_POS"] = df["chr"] + "_" + df["TF_start"].map(str) + "_" + df["TF_end"].map(str)
    df["MOTIF_LEN"] = df["TF_end"] - df["TF_start"] + 1
    df["TF_Name"] = df["MOTIF_ID"].map(y)

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

def sort_and_merge_reference(x):
    for i in x:
        print ("Sorting file: " + i)

        df_sort = pd.read_csv(i, header=None, usecols=[0, 1, 2], low_memory=False)

        df_sort.sort_values(by=[0, 1, 2], inplace=True)

        a = pybedtools.BedTool.from_dataframe(df_sort)
        c = a.merge()
        d = c.moveto(i)

    return x
