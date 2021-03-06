#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import argparse
import pandas as pd
import pybedtools
import math
import numpy as np
import scipy
import glob
from itertools import groupby, chain
from collections import OrderedDict
from os import walk
import itertools
import re

def hg19_UCSC_to_Gene_Name(x):


    df_kgXref = pd.read_csv(x, sep="\t", header=0,dtype={'rfamAcc': str, 'tRnaName': str})
    df_UCSC_Gene_Name = pd.DataFrame()
    df_UCSC_Gene_Name["UCSC_ID"] = df_kgXref["#kgID"]
    df_UCSC_Gene_Name["Gene_Symbol"] = df_kgXref["geneSymbol"]
    df_UCSC_Gene_Name.set_index("UCSC_ID")

    dict_UCSC_Gene_Name = dict(zip(df_UCSC_Gene_Name.UCSC_ID, df_UCSC_Gene_Name.Gene_Symbol))

    return dict_UCSC_Gene_Name

def mm10_UCSC_to_Gene_Name(x):
    df_mm10_UCSC_Gene_Name = pd.read_csv(x, sep = "\t", header = 0)
    columns = ["Gene_Name", "UCSC_Stable_ID"]
    df_mm10_UCSC_Gene_Name.columns = columns

    new_df = pd.DataFrame()
    new_df["UCSC_Stable_ID"] = df_mm10_UCSC_Gene_Name["UCSC_Stable_ID"]
    new_df["Gene_Name"] = df_mm10_UCSC_Gene_Name["Gene_Name"]
    new_df.set_index("UCSC_Stable_ID")

    dict_mm10_GSM_Gene_Name = dict(zip( new_df.UCSC_Stable_ID, new_df.Gene_Name))
    return dict_mm10_GSM_Gene_Name


def GEO_to_TF(y):
    #read the CHIP seq data file as csv
    #this is the file given by Sreeja
    df_GEO = pd.read_csv(y, sep="\t", header=0)
    #Make column names and add them to the df
    columns = ["GSE_ID", "GSM_ID", "SRR#", "CELL_TYPE", "MOLECULE"]
    df_GEO.columns = columns

    #Add new column to df_GEO, splittling the MOLECULE column by : and picking the last string element
    #this give us the TF
    df_GEO["TF"] = df_GEO["MOLECULE"].str.split(":").str[3]

    #Use code similar to this to merge data on known TFs that we are interested in
    #common_GSM_id = pd.merge(df_GEO, df_GSM_id_filename, how='inner', on=['GSM_ID'])

    ##Pull only the GSM_IDs and TFs columns
    ##df_GSM_ID_TF = common_GSM_id.loc[:, ['GSM_ID', 'TF']]

    #pull the GSM_ID and TF into a new df
    df_GSM_ID_TF = pd.DataFrame()
    df_GSM_ID_TF["GSM_ID"] = df_GEO["GSM_ID"]
    df_GSM_ID_TF["TF"] = df_GEO["TF"]
    df_GSM_ID_TF.set_index("GSM_ID")

    #Known tf given by emily and this is hardcoded
    Known_TF = pd.read_csv('/Users/war9qi/Faiz_Workspace/meta/tf_geneSymbol_human_weirauch.txt')
    columns = ["TF"]
    Known_TF.columns = columns
    #Merge dfs by same TF
    df_GSM_ID_TF = pd.merge(df_GSM_ID_TF, Known_TF, how = 'inner')

    dict_GSM_TF = dict(zip(df_GSM_ID_TF.GSM_ID, df_GSM_ID_TF.TF))

    return dict_GSM_TF


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

def bin_parse(x, list_of_bins):
    peaks = x["ID"]
    bins = x['bins']
    start = x['Start']
    chrom = str(x['Chr'])
    counter = 0

    r = pd.Series(range(int(bins)))

    r.apply(bin_write, args=(list_of_bins, chrom, peaks, bins, start))

def bin_write(x, list_of_bins, chrom, peaks, bins, start):
    bin_track = x * 100
    bin_start = start + bin_track
    bin_end = bin_start + 99

    d = pd.DataFrame([chrom, bin_start, bin_end, peaks, bins]).transpose()
            
    list_of_bins.append(d)

def genemeta2dict(x):
    # Create the dataframe that will be used to get the dictionary of Gene names
    df_gene = pd.read_csv(x, sep="\t", header=0)

    # Get the base names and set as the index ---- This might need to be changed later
    df_gene_names = pd.DataFrame()
    df_gene_names["UCSC"] = df_gene["#mm10.knownToEnsembl.name"]
    df_gene_names["Symbol"] = df_gene["mm10.kgXref.geneSymbol"]
    df_gene_names.set_index("UCSC")

    dict_gene = dict(zip(df_gene_names.UCSC, df_gene_names.Symbol))

    return dict_gene

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

def merge_rep_BEDS(x, d, og_dir, modes):
    df = d[d.TF_Name == x]

    for m in modes:

        t_df = df[df.MODE == m]

        tf_file_name = x + "_" + m + ".bed"

        allFiles = t_df["file_path"]

        frame = pd.DataFrame()
        list_ = []

        for ind_file in allFiles:
            con_df = pd.read_csv(ind_file, sep="\t", header=None, usecols=[0,1,2,3], low_memory = False)
            list_.append(con_df)

        os.chdir(og_dir)

        frame = pd.concat(list_)
        frame["TF_name"] = x
        frame.sort_values(by=[0, 1, 2], inplace=True)
        
        a = pybedtools.BedTool.from_dataframe(frame)
        c = a.merge()
        e = c.to_dataframe()
        e["TF_name"] = x
        e.to_csv(tf_file_name, sep="\t", index=False, header=False)

def TFmeta2dict(x):
    # Create the dataframe that will be used to get the dictionary of TF names
    df_TF = pd.read_csv(x, sep="\t", header=0)

    # Get the base names and set as the index ---- This might need to be changed later
    df_TF_names = pd.DataFrame()
    df_TF_names["Motif"] = df_TF["Motif_ID"].str.replace("_1.97d", "")
    df_TF_names["TF"] = df_TF["TF_Name"]
    df_TF_names.set_index("Motif")

    dict_TF = dict(zip(df_TF_names.Motif, df_TF_names.TF))
    
    return dict_TF

def parse_moods(x, y):
    columns = ['ID', 'Motif', 'TF_pos', 'Motif_sequence']
    col_types = {"ID": "object", "Motif": "object", "TF_pos": "uint16", "Motif_sequence": "object"}

    df_moods = pd.read_csv(x, sep="|", usecols=[0,1,2,5], header=None, names=columns, dtype=col_types, low_memory=False)

    df_moods["Motif"] = df_moods["Motif"].str.split('_').str[0]

    df_moods["TF_Name"] = df_moods["Motif"].map(y)

    df_moods["Chr"], df_moods["Coord"] = df_moods["ID"].str.split(':').str

    #df_moods["Start"], df_moods["Stop"] = df_moods["Coord"].str.split("-").str

    #df_moods["TF_start"] = df_moods["Start"].apply(int) + 1 + df_moods["TF_pos"].apply(int)
    df_moods["TF_start"] = df_moods["Coord"].str.split("-").str[0].apply(int) + 1 + df_moods["Coord"].str.split("-").str[1].apply(int)

    df_moods["TF_end"] = df_moods["TF_start"] + df_moods["Motif_sequence"].str.len() - 1

    order= ["Chr", "TF_start", "TF_end", "TF_Name", "ID"]

    df_moods = df_moods[order]

    bed_moods = pybedtools.BedTool.from_dataframe(df_moods)
    
    return bed_moods

def generate_prior(BED, GROUPER, FILENAME, OUT_DIR, FILENAME_BINARY):
    basename = os.path.splitext(FILENAME)[0]
    out_sparse = OUT_DIR + "/" + basename + "_sp.txt"
    out_merged = OUT_DIR + "/" + basename + '_merged.tsv'
    out_merged_sparse = OUT_DIR + "/" + basename + '_merged_sp.tsv'
    merged_filename = OUT_DIR + "/" + basename + '_merged.meta'
    out_sparse_bin = OUT_DIR + "/" + basename + "_binary_sp.txt"

    df_gb = BED.groupby(GROUPER).count()
    df_gb.reset_index(inplace=True)

    df_piv = df_gb.pivot_table(index="Gene_Symbol", columns="TF_name", values="Gene_name", aggfunc=np.sum)
    df_piv.fillna(0, inplace=True)
    
    #print(df_piv)
    
    #df_bin = df_piv.apply(lambda x: [y if y <= 1 else 1 for y in x])
    #print(df_bin)
    #df_bin.to_csv(FILENAME_BINARY, sep="\t")

    df_piv.to_csv(FILENAME, sep="\t")
    
    df_prior = df_piv.unstack().reset_index(name="Weight")
    df_prior.columns = ["regulator", "target", "weight"]
    df_prior = df_prior[df_prior.weight != 0]
    df_prior.to_csv(out_sparse, sep="\t", index=False)
    
    #df_prior_bin = df_bin.unstack().reset_index(name="Weight")
    #df_prior_bin.columns = ["regulator", "target", "weight"]
    #df_prior_bin = df_prior_bin[df_prior_bin.weight != 0]
    #df_prior_bin.to_csv(out_sparse_bin, sep="\t", index=False)
    
    print ("################## sum ##################")
    print (str(len(df_prior.index)) + ' interactions.')
    print (str(len(df_prior["target"].unique())) + ' genes with at least one TF interaction.')
    print (str(len(df_prior["regulator"].unique())) + ' total TFs in prior.')
    print (FILENAME + ' generated.')
    print (out_sparse + ' generated.')
    
    df_trans = df_piv.transpose()
    df_trans["sum"] = df_trans.sum(axis=1)
    df_trans["count"] = 0

    gp_num = df_trans.groupby("sum").count()

    df_merge = gp_num[gp_num["count"] > 1]
    df_merge.reset_index(inplace=True)

    keep = list(df_merge["sum"])

    df_merge_names = df_trans.loc[df_trans['sum'].isin(keep)]
    df_merge_names["is_dup"] = df_merge_names.duplicated(keep=False)
    df_merge_names = df_merge_names[df_merge_names["is_dup"] == True]
    
    merge_list = df_merge_names["sum"].unique().tolist()
    
    if len(merge_list) >= 1:
        for i in merge_list:
            merge_df = df_merge_names.loc[df_merge_names['sum'] == i]
            merge_list_tfs = merge_df.index.tolist()
            
            if len(merge_list_tfs) >= 3:
                final_name = merge_list_tfs[0] + "_" + merge_list_tfs[1] + "..."
                
            else:
                final_name = merge_df.index.str.cat(sep="_")
            
            df_piv.rename(index=str, columns={merge_list_tfs[0]: final_name}, inplace=True)
            df_piv.drop(labels=merge_list_tfs[1:], axis=1, inplace=True)
    
    else:
        pass

    df_prior = df_piv.unstack().reset_index(name="Weight")
    df_prior.columns = ["regulator", "target", "weight"]
    df_prior = df_prior[df_prior.weight != 0]
    
    df_piv.to_csv(out_merged, sep="\t")

    df_prior.to_csv(out_merged_sparse, sep="\t", index=False)
   
    print ("################## sum - merged ##################")
    print (str(len(df_prior.index)) + ' interactions.')
    print (str(len(df_prior["target"].unique())) + ' genes with at least one TF interaction.')
    print (str(len(df_prior["regulator"].unique())) + ' total TFs in prior.')
    print (out_merged + ' generated.')
    print (out_merged_sparse + ' generated.')

def pybedtools_intersect(x, y):
    bed_1 = pybedtools.BedTool(x)
    bed_2 = pybedtools.BedTool(y)

    bed_intersect = bed_1.intersect(bed_2, wa=True, wb=True)

    return bed_intersect

def priorTable2Sparse(x, y, z):
    """x = file, y=basename, z=output_dir"""
    out_sparse = z + "/" + y + "_sp.txt"
    
    df = pd.read_csv(x, sep="\t", header=0, index_col=0)

    df_prior = df.unstack().reset_index(name="Weight")

    df_prior.columns = ["regulator", "target", "weight"]
    
    df_prior = df_prior[df_prior.weight != 0]

    gb_TF = df_prior.groupby("regulator").count()

    gb_TF.drop("weight", axis=1, inplace=True)

    gb_TF.reset_index(inplace=True)

    df_prior.to_csv(out_sparse, sep="\t", index=False)

    print (str(len(df_prior.index)) + ' interactions.')
    print (str(len(df_prior["target"].unique())) + ' genes with at least one TF interaction.')
    print (str(len(df_prior["regulator"].unique())) + ' total TFs in prior.')
    print (out_sparse + ' generated.')

    return df

def mergeDegeneratePriorTFs(x, y):
    out_merged_sparse = y + '_merged_sp.tsv'
    out_merged = y + '_merged.tsv'

    df_trans = x.transpose()

    df_trans["sum"] = df_trans.sum(axis=1)
    
    df_trans["count"] = 0

    gp_num = df_trans.groupby("sum").count()

    df_merge = gp_num[gp_num["count"] > 1]
    
    df_merge.reset_index(inplace=True)

    keep = list(df_merge["sum"])

    df_merge_names = df_trans.loc[df_trans['sum'].isin(keep)]

    df_merge_names.sort_values("sum")

    df_merge_names["is_dup"] = df_merge_names.duplicated(keep=False)
    
    tf_to_merge = df_merge_names[df_merge_names["is_dup"] == True]

    if len(tf_to_merge) == 0:
        print("No TFs to merge")
    
    else:
        x.to_csv(out_merged_sparse, sep="\t")

class BINS:
    def __init__(self, path, DHS_BED, bin_dir, ofd, ChIP_dir):
        self.DHS_BED = os.path.abspath(DHS_BED)
        self.ofd = ofd
        
        f_ = []

        for root,dirs,filenames in os.walk(path):
            filenames = [g for g in filenames if not g[0] == '.']
            dirs[:] = [d for d in dirs if not d[0] == '.']
            for f in filenames:
                relFile = os.path.abspath(os.path.join(root, f))
                f_.append(relFile)

        f_.append(self.DHS_BED)
        self.files = f_

        self.bin_dir = bin_dir

    def bin_group_collect(self, ofd, bin_DHS, blacklist):        
        if bin_DHS == True:
            f_name = "BINNED_DHS_ONLY.bed"
            final_name = "DHS.bin.bed"
            DHS_BED = self.DHS_BED
            self.dhs = final_name

            a = pybedtools.BedTool(DHS_BED)
            b = pybedtools.BedTool(blacklist)

            a_and_b = a.intersect(b, v=True)     
            
            df = a_and_b.to_dataframe(names=["Chr", "Start", "Stop"])

        else:
            f_name = "BINNED_UNION.bed"
            final_name = "Union.bin.bed"
            list_df = []
            self.union = final_name

            for fname in self.files:
                bin_df = pd.read_csv(fname, sep="\t", header=None, usecols=[0,1,2])
                list_df.append(bin_df)
            
            frame = pd.concat(list_df)

            frame.sort_values(by=[0, 1, 2], inplace=True)

            a = pybedtools.BedTool.from_dataframe(frame)
            c = a.merge()
        
            df = c.to_dataframe(names=["Chr", "Start", "Stop"])
        
        df["len"] = df["Stop"] - df["Start"] + 1
        df["bins"] = df["len"]/100
        df["bins"] = df["bins"].apply(math.ceil)
        df["ID"] = df["Chr"] + "_" + df["Start"].apply(str) + "_" + df["Stop"].apply(str)
        
        os.chdir(self.bin_dir)

        df.to_csv(f_name, sep="\t", index=False)

        list_of_bins = []

        df.apply(bin_parse, axis=1, args=(list_of_bins,)) 
        
        df = pd.concat(list_of_bins)
        df['index']=df.reset_index().index
        df["bin_ID"] = "bin_" + df["index"].apply(str)
        df.to_csv(final_name, sep="\t", index=False, header=False)
        
        os.chdir(self.ofd)

    def intersect_chip_bin(self, modes, bin_file, percentile_folder, intersect_dir, bin_name, tf_names):
        print ("Intersect ChIP with BINS")
        for i in percentile_folder:
            base = os.path.basename(i)
            
            for m in modes:
                os.chdir(i)
                list_df = []
                
                for fname in glob.glob("*" + m + ".bed"):
                    bin_df = pd.read_csv(fname, sep="\t", header=None)
                    list_df.append(bin_df)
                
                frame = pd.concat(list_df)

                a = pybedtools.BedTool((self.bin_dir + "/" + bin_file))
                b = pybedtools.BedTool.from_dataframe(frame)

                a_and_b = a.intersect(b, wa=True, wb=True)

                os.chdir(intersect_dir)
                make_set_dir(bin_name, False)

                f_name = bin_name + "_" + base + "_" + m + ".bed"
                c = a_and_b.moveto(f_name)

    def intersect_moods_bin(self, bin_file, MOODS_dir, intersect_dir, n, tf_names):
        os.chdir(MOODS_dir)
        y = glob.glob("*.bed")
        for i in y:
            p_val = i.replace("GM12878_DHS_MOTIFS_", "").replace(".bed", "")

            os.chdir(MOODS_dir)

            a = pybedtools.BedTool((self.bin_dir + "/" + bin_file))
            b = pybedtools.BedTool(i)
            a_and_b = a.intersect(b, wa=True, wb=True)

            os.chdir(intersect_dir)
            make_set_dir(n, False)

            c = a_and_b.moveto(n + "_MOODS_" + p_val + ".bed")
            
class MARIO:
    def __init__(self, chip_bed, percentile, blacklist):
        self.chip_bed = chip_bed
        self.percentile = percentile
        self.blacklist = blacklist

    def parse_singles_percentile(self, df, percentile):
        df["file_path"].apply(parse_percentile, args = (df, percentile, True, self.blacklist))

    def parse_replicate_percentile(self, df, percentile, og_dir, modes, tf_df):
        df["file_path"].apply(parse_percentile, args = (df, percentile, False, self.blacklist))

        tf = pd.DataFrame(tf_df)
        tf[0].apply(merge_rep_BEDS, args = (df, og_dir, modes))

class META:
    def __init__(self, META, chip_bed, tf_name_file):
        self.META_data = META
        self.chip_bed = chip_bed
        self.exclude_list = ["viral_protein", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3", "H3K79me2"]
    
        # Create an empty dataframe
        df_dict = pd.DataFrame()

        # Create a DF from TF_FILE_associations
        TF_df = pd.read_csv(tf_name_file, sep="\t", header=0)
        df_dict["Motif"] = TF_df["Motif_ID"]
        df_dict["TF"] = TF_df["TF_Name"]
        df_dict.set_index("Motif")

        #Create a dictionary of all TF and Motif associations
        dicted = dict(zip(df_dict.Motif, df_dict.TF))

        self.tf_dict = dicted

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

        unique_TFS = CHIP_df["TF_Name"].unique()

        self.unique_TF_single_names = unique_TF_single
        self.unique_TF_reps_names = unique_TF_reps
        self.unique_MODES = unique_MODES

        self.unique_tf_single_df = unique_tf_single_df
        self.unique_tf_rep_df = unique_tf_rep_df
        self.chip_df = CHIP_df

        self.unique_TFS = unique_TFS

class MOODS:
    """This object is being used to store the information and process
    the data data is associated with each MOODs file at various cutoffs
    This will take as inputs:
        1) file_path
        2) tf_name_file
        3) output_dir"""

    # Initialize class with these parameters
    def __init__(self, file_path, output_dir, domain_bed, tf_dict, unique_TFS):
        self.file_path = file_path
        self.output_file_directory = output_dir
        self.column_names = ['PEAK_ID', 'PWM_FILE', 'TF_START', 'MOTIF_SEQ']
        self.column_types = {
            'PEAK_ID':          "category",
            'PWM_FILE':         "category",
            'TF_START':         "uint16",
            'MOTIF_SEQ':        "category"}
       
        self.domain_bed = domain_bed

        # Create a DF of the MOODs file, setting dtypes, columns, and names
        df = pd.read_csv(self.file_path, usecols=[0,1,2,5], names=self.column_names, dtype=self.column_types, header=None, sep="|", low_memory=False)

        #Parse through the file and extract names, frames, etc
        df['PWM_FILE'] = df['PWM_FILE'].str.replace('_JASPAR.txt.pfm', '')
        
        df["TF_Name"] = df["PWM_FILE"].map(tf_dict)

        df = df[df['TF_Name'].isin(unique_TFS)]

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

        #Not sure if the below code is faster or if copy the DF over itself is faster
        #df.drop(columns=["TF_START", "PWM_FILE", "MATCH_SCORE", "MOTIF_POS", "MOTIF_LEN", "MOTIF_SEQ", "STRAND", "start", "stop"], inplace=True)

        df = df[["chr", "TF_start", "TF_end", "PEAK_ID", "TF_Name"]]

        MOODS_HITS_BED = os.path.basename(self.file_path)
        MOODS_HITS_BED = MOODS_HITS_BED.replace(".txt", ".bed")

        df.sort_values(by=['chr', "TF_start", "TF_end"], inplace=True)

        self.moods_bed = MOODS_HITS_BED

        a = pybedtools.BedTool.from_dataframe(df)

        del df

        b = pybedtools.BedTool(self.domain_bed)

        a_and_b = a.intersect(b, wa=True, wb=True)

        c = a_and_b.moveto(self.moods_bed)

        df = c.to_dataframe()

        df["SUB_ID"] = df[5] + "_" + df[6].apply(str) + "_" + df[7].apply(str)

        MOODS_HITS_FN = self.moods_bed.replace(".txt", ".bed")

        df.to_csv(MOODS_HITS_FN, sep="\t", header = False, index=False, columns= [0, 1, 2, 4, 3, 8, "SUB_ID"])

        self.moods_bed_intersect =  MOODS_HITS_FN

class RESULTS:
    def __init__(self):
        self.column_types = {
            'BIN_ID': "category",
            'TF_NAME': "category"}

    def parse_MOODS(self):
        MOODS_file_list = glob.glob("*.bed")
        MOODS_df_list = []

        for i in MOODS_file_list:
            df = pd.read_csv(i, sep="\t", 
                            header=None,
                            usecols=[6, 10],
                            dtype=self.column_types,
                            names=self.column_types.keys())
            df["SET_ID"] = i
            MOODS_df_list.append(df)

        MOODS_con = pd.concat(MOODS_df_list)

        del MOODS_df_list

        gbm = MOODS_con.groupby(["SET_ID", df.BIN_ID.astype(object), df.TF_NAME.astype(object)]).count()

        return gbm

    def parse_CHIP(self):
        CHIP_FILE_LIST = glob.glob("*.bed")
        CHIP_DF_LIST = []

        for i in CHIP_FILE_LIST:
            df = pd.read_csv(i, sep="\t", 
                            header=None,
                            usecols=[6, 10],
                            dtype=self.column_types,
                            names=self.column_types.keys())
            df["SET_ID"] = i
            CHIP_DF_LIST.append(df)

        CHIP_con = pd.concat(CHIP_DF_LIST)
        
        del CHIP_DF_LIST

        gbc = CHIP_con.groupby(["SET_ID", df.BIN_ID.astype(object), df.TF_NAME.astype(object)]).count()

        return gbc
