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
from os import walk
import time

#COLOR PARAMETERS
moods_color = colored.fg(202) + colored.attr(1)
MARIO_PARSE_color = colored.fg(14) + colored.attr(1)
MAKE_DIR_color = colored.fg(201) + colored.attr(1)
tacman_color = colored.fg(226) + colored.attr(1)
checkpoint = stylize("################################################################################################################", tacman_color)

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

def clock(x, y, z):
    elapsed = (time.time() - x)
    if elapsed <= 60:
            print (stylize("---It took TACMAN %s seconds---" % elapsed, tacman_color))

    else:
            minutes = elapsed/60
            print (stylize("---It took TACMAN %s minutes---" % minutes, tacman_color))

    elapsed_min = (time.time() - y)/60
    print (stylize("---TACMAN has been running for %s minutes ---" % elapsed_min, tacman_color))

def make_set_dir(x, y):
    cwd = os.getcwd()

    dir = cwd + "/" + x

    if not os.path.exists(dir):
            os.makedirs(dir)

    os.chdir(dir)

    if y == True:
        return dir
    else:
        pass

def merge_rep_BEDS(x, d, og_dir, modes):
    df = d[d.TF_Name == x]
    TF_n = d["TF_Name"].unique()

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
        frame["TF_name"] = TF_n
        frame.sort_values(by=[0, 1, 2], inplace=True)
        frame.to_csv(tf_file_name, sep="\t", index=False, header=False)

        a = pybedtools.BedTool(tf_file_name)
        c = a.merge()
        d = c.moveto(tf_file_name)

def parse_percentile(x, y, p, tf_name, blacklist):
    TF_name = y.loc[y['file_path'] == x, 'TF_Name'].iloc[0]
    MODE_name = y.loc[y['file_path'] == x, 'MODE'].iloc[0]
    og_file_name = os.path.basename(x)
    tf_file_name = TF_name + "_" + MODE_name + ".bed"

    a = pybedtools.BedTool(x)
    b = pybedtools.BedTool(blacklist)
    a_and_b = a.intersect(b, v=True)     
    df = a_and_b.to_dataframe()
    
    col_num = len(df.columns)
    
    if col_num == 11:
        per_len = len(df.index)/p
        df.sort_values(by=["thickEnd"], ascending=False, inplace=True)
        df = df.head(per_len)
        df["TF_name"] = TF_name

        if tf_name == True:
            df.to_csv(tf_file_name, sep="\t", index=False, header = False, columns= ["chrom","start","end","TF_name"])

        else:
            df.to_csv(og_file_name, sep="\t", index=False, header = False, columns= ["chrom","start","end","TF_name"])

    elif col_num == 13:
        per_len = len(df.index)/p
        df.sort_values(by=[12], ascending=False, inplace=True)
        df = df.head(per_len)
        df["TF_name"] = TF_name

        if tf_name == True:
            df.to_csv(tf_file_name, sep="\t", index=False, header = False, columns= [0,1,2,"TF_name"])

        else:
            df.to_csv(og_file_name, sep="\t", index=False, header = False, columns= [0,1,2,"TF_name"])

class BINS:
    def __init__(self, path, DHS_list, bin_dir):
        f = []
        
        for dir_, _, files in os.walk(path):
            for fileName in files:
                relDir = os.path.relpath(dir_, path)
                relFile = os.path.join(relDir, fileName)
                f.append(relFile)
        
        self.ofd = path
        self.files = f
        self.DHS_list = DHS_list
        self.bin_dir = bin_dir

    def bin_group_collect(self, ofd, bin_DHS):        
        binned_genome_meta = []

        if bin_DHS == True:
            f_name = "BINNED_DHS_ONLY.bed"
            final_name = "DHS.bin"
            fname = self.DHS_list
            frame = pd.read_csv(fname, sep="\t", header=None, usecols=[0,1,2,3])

        else:
            f_name = "BINNED_UNION.bed"
            final_name = "Union.bin"
            files_to_bin = self.files[1:]
            list_df = []

            for fname in files_to_bin:
                bin_df = pd.read_csv(fname, sep="\t", header=None, usecols=[0,1,2,3])
                list_df.append(bin_df)
            
            frame = pd.concat(list_df)

        frame.sort_values(by=[0, 1, 2], inplace=True)

        os.chdir(self.bin_dir)

        a = pybedtools.BedTool.from_dataframe(frame)
        c = a.merge()
        df = c.moveto(f_name)
        
        df = pd.read_csv(f_name, sep="\t", header=None, usecols=[0,1,2,3])

        df.columns = ["Chr", "Start", "Stop"]

        df["len"] = df["Stop"] - df["Start"] + 1
        df["bins"] = df["len"]/100
        df["bins"] = df["bins"].apply(math.ceil)
        df["ID"] = df["Chr"] + "_" + df["Start"].apply(str) + "_" + df["Stop"].apply(str)

        df.to_csv(f_name, sep="\t", index=False)

        binned_file = df.apply(bin_parse, axis = 1) 

        binned_file.to_csv(final_name, sep="\t", index=False)

        os.chdir(self.ofd)

        bin_files = []
        self.bin_files = bin_files
        self.bin_files.append(final_name)

    def intersect_bin(self, modes):
        for m in modes:
            y = glob.glob("*" + m + ".bed")

            for i in self.bin_files:
                print ("Intersect BED: List to intersect: " + y)

                a = pybedtools.BedTool(i)
                b = pybedtools.BedTool(y)

                print ("Intersect BED: Working on: " + f_name)

                a_and_b = a.intersect(b, wa=True, wb=True)

                c = a_and_b.moveto("ChIP_BINS" + m + ".bed")

                d = pybedtools.BedTool(x)

                print ("Intersect BED: Working on: CHIP")

                a_and_d = a.intersect(d, wa=True, wb=True)

                e = a_and_d.moveto("MOODS_DHS_RE_SUB_BINS_" + m + ".bed")

class MARIO:
    def __init__(self, chip_bed, percentile, blacklist):
        self.chip_bed = chip_bed
        self.percentile = percentile
        self.blacklist = blacklist

    def parse_singles_percentile(self, df, percentile):
        percentiles = 100/percentile
        df["file_path"].apply(parse_percentile, args = (df, percentile, True, self.blacklist))

    def parse_replicate_percentile(self, df, percentile):
        percentiles = 100/percentile
        df["file_path"].apply(parse_percentile, args = (df, percentile, False, self.blacklist))

    def merge_replicate_BEDS(self, df, og_dir, modes, tf_df):
        tf = pd.DataFrame(tf_df)
        tf[0].apply(merge_rep_BEDS, args = (df, og_dir, modes))

class META:
    def __init__(self, META, chip_bed):
        self.META_data = META
        self.chip_bed = chip_bed
        self.exclude_list = ["viral_protein", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3", "H3K79me2"]

    def META_parse(self):
        df = pd.read_csv(self.META_data, sep="\t", header=None, usecols=[1, 13], names=["Sample_Name", "info"])

        df["TF"] = df["info"].str.split(":").str[3]
        df["Type"] = df["info"].str.split(":").str[1]

        for i in self.exclude_list:
            df = df[df.Type != i]
            df = df[df.TF != i]

        dicted = dict(zip(df.Sample_Name, df.TF))

        CHIP_df = pd.DataFrame.from_dict(self.chip_bed)
        CHIP_df["Basename"] = CHIP_df[0].apply(os.path.basename)
        CHIP_df["File_ext"] = CHIP_df.Basename.str.split("_").str[-1]
        CHIP_df["MODE"] = CHIP_df.File_ext.str.split(".").str[0]
        CHIP_df["Sample_Name"] = CHIP_df.Basename.str.split("_").str[0]
        CHIP_df["TF_Name"] = CHIP_df["Sample_Name"].map(dicted)
        CHIP_df["file_path"] = CHIP_df[0]
        CHIP_df.dropna(inplace=True)

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
        a = pybedtools.BedTool(self.moods_bed)
        b = pybedtools.BedTool(self.domain_bed)

        a_and_b = a.intersect(b, wa=True, wb=True)

        c = a_and_b.moveto(self.moods_bed)

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

        x = self.intersect_df[[0, 1, 2, 4, 3, 8, "SUB_ID"]]

        MOODS_HITS_FN = self.moods_bed.replace(".txt", ".bed")

        x.to_csv(MOODS_HITS_FN, sep="\t", header = False, index=False)

        self.moods_bed_intersect =  MOODS_HITS_FN

    def parse_moods(self):
        # Create a DF of the MOODs file, setting dtypes, columns, and names
        df = pd.read_csv(self.file_path, usecols=[0,1,2,3,4,5], names=self.column_names, dtype=self.column_types, header=None, sep="|")

        #Parse through the file and extract names, frames, etc
        df['PWM_FILE'] = df['PWM_FILE'].str.replace('_JASPAR.txt.pfm', '')

        df_tmp1 = df['PEAK_ID'].str.split(":", expand=True)
        df_tmp2 = df_tmp1[1].str.split("-", expand=True)

        df.drop(columns=['PEAK_ID'], inplace=True)

        df["chr"] = df_tmp1[0]
        df["start"] = df_tmp2[0]
        df["stop"] = df_tmp2[1]

        #delete frames to save space
        del df_tmp1
        del df_tmp2

        df["TF_start"] = df["start"].apply(int) + 1 + df["TF_START"]
        df["TF_end"] = df["TF_start"] + df["MOTIF_SEQ"].str.len() - 1

        df["PEAK_ID"] = df["chr"] + "_" + df["start"].apply(str) + "_" + df["stop"].apply(str)
        df["MOTIF_POS"] = df["chr"] + "_" + df["TF_start"].apply(str) + "_" + df["TF_end"].apply(str)

        df["MOTIF_LEN"] = df["TF_end"] - df["TF_start"] + 1
        df["TF_Name"] = df["PWM_FILE"].map(self.tf_dict)

        #Not sure if the below code is faster or if copy the DF over itself is faster
        #df.drop(columns=["TF_START", "PWM_FILE", "MATCH_SCORE", "MOTIF_POS", "MOTIF_LEN", "MOTIF_SEQ", "STRAND", "start", "stop"], inplace=True)

        df = df[["chr", "TF_start", "TF_end", "PEAK_ID", "TF_Name"]]

        MOODS_HITS_BED = os.path.basename(self.file_path)
        MOODS_HITS_BED = MOODS_HITS_BED.replace(".txt", ".bed")

        df.sort_values(by=['chr', "TF_start", "TF_end"], inplace=True)

        df.to_csv(MOODS_HITS_BED, sep="\t", index=False, header=False)

        self.moods_bed = MOODS_HITS_BED
