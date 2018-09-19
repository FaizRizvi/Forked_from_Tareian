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

def dict_TF_df(x):
	df = pd.DataFrame()

	df["Motif"] = x["Motif_ID"]
	df["TF"] = x["TF_Name"]

	df.set_index("Motif")

	dicted = dict(zip(df.Motif, df.TF))

	return dicted

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

def bed_intersect(x, y):
    a = pybedtools.BedTool(x)
    b = pybedtools.BedTool(y)

    a_and_b = a.intersect(b, wa=True, wb=True)

    c = a_and_b.moveto("TF_domain.bed")

    df = c.to_dataframe()

    return df
