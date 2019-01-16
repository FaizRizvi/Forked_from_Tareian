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
    
    df_piv.to_csv(FILENAME, sep="\t")
    
    df_prior = df_piv.unstack().reset_index(name="Weight")
    df_prior.columns = ["regulator", "target", "weight"]
    df_prior = df_prior[df_prior.weight != 0]
    df_prior.to_csv(out_sparse, sep="\t", index=False)
    
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
