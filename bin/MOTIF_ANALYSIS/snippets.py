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

#COLOR PARAMETERS
moods_color = colored.fg(202) + colored.attr(1)
MARIO_PARSE_color = colored.fg(14) + colored.attr(1)
MAKE_DIR_color = colored.fg(201) + colored.attr(1)

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

def make_set_dir(x, y):
    cwd = os.getcwd()

    dir = cwd + "/" + x
    print (stylize("Making Directory: " + dir, MAKE_DIR_color))

    if not os.path.exists(dir):
            os.makedirs(dir)

    print (stylize("Making Directory: Changing working directory to: " + dir, MAKE_DIR_color))
    os.chdir(dir)

    if y == True:
        return dir
    else:
        pass

def merge_rep_BEDS(x, d, og_dir, modes):
    print ("Merging Replicate BEDS for TF: " + x)
    df = d[d.TF_Name == x]

    for m in modes:

        t_df = df[df.MODE == m]

        tf_file_name = x + "_" + m + ".bed"

        allFiles = t_df["file_path"]

        frame = pd.DataFrame()
        list_ = []

        for ind_file in allFiles:
            con_df = pd.read_csv(ind_file, sep="\t", header=None, usecols=[0,1,2,3])
            list_.append(con_df)

        os.chdir(og_dir)

        frame = pd.concat(list_)

        frame.sort_values(by=[0, 1, 2], inplace=True)
        frame.to_csv(tf_file_name, sep="\t", index=False, header=False)

        a = pybedtools.BedTool(tf_file_name)
        c = a.merge()
        d = c.moveto(tf_file_name)

def parse_percentile(x, y, p, tf_name):
    TF_name = y.loc[y['file_path'] == x, 'TF_Name'].iloc[0]
    MODE_name = y.loc[y['file_path'] == x, 'MODE'].iloc[0]
    og_file_name = os.path.basename(x)
    tf_file_name = TF_name + "_" + MODE_name + ".bed"

    df = pd.read_csv(x, sep="\t", header=None)
    col_num = len(df.columns)

    if col_num == 11:
        per_len = (len(df[7]))/p
        df.sort_values(by=[7], ascending=False, inplace=True)
        df = df.head(per_len)

        if tf_name == True:
            df.to_csv(tf_file_name, sep="\t", index=False, header = False, columns= [0,1,2,3])

        else:
            df.to_csv(og_file_name, sep="\t", index=False, header = False, columns= [0,1,2,3])

    elif col_num == 13:
        per_len = (len(df[12]))/p
        df.sort_values(by=[12], ascending=False, inplace=True)
        df = df.head(per_len)

        if tf_name == True:
            df.to_csv(tf_file_name, sep="\t", index=False, header = False, columns= [0,1,2,3])

        else:
            df.to_csv(og_file_name, sep="\t", index=False, header = False, columns= [0,1,2,3])

class BINS:
    def __init__(self, modes, MOODS):
        self.modes = modes
        self.MOODS = MOODS

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

class MARIO:
    def __init__(self, chip_bed, percentile):
        self.chip_bed = chip_bed
        self.percentile = percentile

    def parse_singles_percentile(self, df, percentile):
        print ("Parsing the " + str(percentile) + " percentile from single files")
        df["file_path"].apply(parse_percentile, args = (df, percentile, True))

    def parse_replicate_percentile(self, df, percentile):
        print ("Parsing the " + str(percentile) + " percentile from replicate files")
        df["file_path"].apply(parse_percentile, args = (df, percentile, False))

    def merge_replicate_BEDS(self, df, og_dir, modes, tf_df):
        print ("Merging replicate files in " + og_dir)
        tf = pd.DataFrame(tf_df)

        tf[0].apply(merge_rep_BEDS, args = (df, og_dir, modes))

class META:
    def __init__(self, META, chip_bed):
        self.META_data = META
        self.chip_bed = chip_bed
        self.exclude_list = ["viral_protein", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3", "H3K79me2"]

    def META_parse(self):
        print (stylize("MARIO Parsing: Parse META Meta: Reading in META Meta File", MARIO_PARSE_color))
        df = pd.read_csv(self.META_data, sep="\t", header=None, usecols=[1, 13], names=["Sample_Name", "info"])

        print (stylize("MARIO Parsing: Parse META Meta: Dropping non-TF ChIP Files", MARIO_PARSE_color))
        df["TF"] = df["info"].str.split(":").str[3]
        df["Type"] = df["info"].str.split(":").str[1]

        for i in self.exclude_list:
            df = df[df.Type != i]
            df = df[df.TF != i]

        print (stylize("MARIO Parsing: Parse META Meta: Mapping Names", MARIO_PARSE_color))
        dicted = dict(zip(df.Sample_Name, df.TF))

        print (stylize("MARIO Parsing: Parse META Meta: Creating ChIP DF name associations based on META Meta File", MARIO_PARSE_color))
        CHIP_df = pd.DataFrame.from_dict(self.chip_bed)
        CHIP_df["Basename"] = CHIP_df[0].apply(os.path.basename)
        CHIP_df["File_ext"] = CHIP_df.Basename.str.split("_").str[-1]
        CHIP_df["MODE"] = CHIP_df.File_ext.str.split(".").str[0]
        CHIP_df["Sample_Name"] = CHIP_df.Basename.str.split("_").str[0]
        CHIP_df["TF_Name"] = CHIP_df["Sample_Name"].map(dicted)
        CHIP_df["file_path"] = CHIP_df[0]
        CHIP_df.dropna(inplace=True)

        print (stylize("MARIO Parsing: Parse META Meta: Creating reference frames", MARIO_PARSE_color))

        #take the input DF and turn it into a groupby object with the count number for each tf and mode
        TF_file_counts = CHIP_df.groupby(["TF_Name", "MODE"]).count()
        TF_file_counts.reset_index(inplace=True)

        #Create a df from the with the sample names that are greater than 1 meaning they need to be merged
        TF_replicates = TF_file_counts[TF_file_counts.Sample_Name > 1]
        TF_singles = TF_file_counts[TF_file_counts.Sample_Name == 1]

        #Parse the percentile from the files and then save them with TF names
        unique_TF_single = TF_singles["TF_Name"].unique()
        unique_TF_reps = TF_replicates["TF_Name"].unique()

        unique_tf_single_df = CHIP_df.loc[CHIP_df['TF_Name'].isin(unique_TF_single)]
        unique_tf_rep_df = CHIP_df.loc[CHIP_df['TF_Name'].isin(unique_TF_reps)]

        unique_MODES = CHIP_df["MODE"].unique()

        self.unique_TF_single_names = unique_TF_single
        self.unique_TF_reps_names = unique_TF_reps
        self.unique_MODES = unique_MODES

        self.unique_tf_single_df = unique_tf_single_df
        self.unique_tf_rep_df = unique_tf_rep_df
        self.chip_df = CHIP_df

class MOODS:
    """This object is being used to store the information and process
    the data data is associated with each MOODs file at various cutoffs
    This will take as inputs:
        1) file_path
        2) tf_name_file
        3) output_dir"""

    # Initialize class with these parameters
    def __init__(self, file_path, tf_name_file, output_dir, domain_bed):
        self.file_path = file_path
        self.tf_name_file = tf_name_file
        self.output_file_directory = output_dir
        self.column_names = ['PEAK_ID', 'PWM_FILE', 'TF_START', 'STRAND', 'MATCH_SCORE', 'MOTIF_SEQ']
        self.column_types = {
            'PEAK_ID': "category",
            'PWM_FILE': "category",
            'TF_START': "uint16",
            'STRAND': "category",
            'MATCH_SCORE': "float32",
            'MOTIF_SEQ': "category"
        }
        self.domain_bed = domain_bed

    def bed_intersect(self):
        print (stylize("Parsing MOODS: Bed Intersect with Domains: Importing BEDS", moods_color))
        a = pybedtools.BedTool(self.moods_bed)
        b = pybedtools.BedTool(self.domain_bed)

        print (stylize("Parsing MOODS: Bed Intersect with Domains: Intersecting BEDS", moods_color))
        a_and_b = a.intersect(b, wa=True, wb=True)

        print (stylize("Parsing MOODS: Bed Intersect with Domains: Writing BEDS", moods_color))
        c = a_and_b.moveto(self.moods_bed)

        print (stylize("Parsing MOODS: Bed Intersect with Domains: Converting to frame", moods_color))
        df = c.to_dataframe()

        self.intersect_df = df

    def dict_TF_df(self):
        # Create an empty dataframe
        df = pd.DataFrame()

        # Create a DF from TF_FILE_associations
        TF_df = pd.read_csv(self.tf_name_file, sep="\t", header=0)
        df["Motif"] = TF_df["Motif_ID"]
        df["TF"] = TF_df["TF_Name"]
        df.set_index("Motif")

        #Create a dictionary of all TF and Motif associations
        dicted = dict(zip(df.Motif, df.TF))

        self.tf_dict = dicted

    def domain_parse(self):
            self.intersect_df["SUB_ID"] = self.intersect_df[5] + "_" + self.intersect_df[6].apply(str) + "_" + self.intersect_df[7].apply(str)

            print (stylize("Parsing MOODS: Ordering frame", moods_color))
            x = self.intersect_df[[0, 1, 2, 4, 3, 8, "SUB_ID"]]

            MOODS_HITS_FN = self.moods_bed.replace(".txt", ".bed")

            print (stylize("Parsing MOODS: Writing to CSV", moods_color))
            x.to_csv(MOODS_HITS_FN, sep="\t", header = False, index=False)

            self.moods_bed_intersect =  MOODS_HITS_FN

    def parse_moods(self):
        # Create a DF of the MOODs file, setting dtypes, columns, and names
        print (stylize("Parsing MOODS: Reading in CSV", moods_color))
        df = pd.read_csv(self.file_path, usecols=[0,1,2,3,4,5], names=self.column_names, dtype=self.column_types, header=None, sep="|")

        #Parse through the file and extract names, frames, etc
        df['PWM_FILE'] = df['PWM_FILE'].str.replace('_JASPAR.txt.pfm', '')

        print (stylize("Parsing MOODS: Splitting and expanding frames", moods_color))
        df_tmp1 = df['PEAK_ID'].str.split(":", expand=True)
        df_tmp2 = df_tmp1[1].str.split("-", expand=True)

        df.drop(columns=['PEAK_ID'], inplace=True)

        df["chr"] = df_tmp1[0]
        df["start"] = df_tmp2[0]
        df["stop"] = df_tmp2[1]

        #delete frames to save space
        del df_tmp1
        del df_tmp2

        print (stylize("Parsing MOODS: Calculating TF motif location", moods_color))
        df["TF_start"] = df["start"].apply(int) + 1 + df["TF_START"]
        df["TF_end"] = df["TF_start"] + df["MOTIF_SEQ"].str.len() - 1

        print (stylize("Parsing MOODS: Extracting Peak and Motif information", moods_color))
        df["PEAK_ID"] = df["chr"] + "_" + df["start"].apply(str) + "_" + df["stop"].apply(str)
        df["MOTIF_POS"] = df["chr"] + "_" + df["TF_start"].apply(str) + "_" + df["TF_end"].apply(str)

        print (stylize("Parsing MOODS: Building Dataframe", moods_color))
        df["MOTIF_LEN"] = df["TF_end"] - df["TF_start"] + 1
        df["TF_Name"] = df["PWM_FILE"].map(self.tf_dict)

        #Not sure if the below code is faster or if copy the DF over itself is faster
        #df.drop(columns=["TF_START", "PWM_FILE", "MATCH_SCORE", "MOTIF_POS", "MOTIF_LEN", "MOTIF_SEQ", "STRAND", "start", "stop"], inplace=True)

        df = df[["chr", "TF_start", "TF_end", "PEAK_ID", "TF_Name"]]

        MOODS_HITS_BED = os.path.basename(self.file_path)
        MOODS_HITS_BED = MOODS_HITS_BED.replace(".txt", ".bed")

        print (stylize("Parsing MOODS: Sorting columns", moods_color))
        df.sort_values(by=['chr', "TF_start", "TF_end"], inplace=True)

        print (stylize("Parsing MOODS: Writing CSV", moods_color))
        df.to_csv(MOODS_HITS_BED, sep="\t", index=False, header=False)

        self.moods_bed = MOODS_HITS_BED
