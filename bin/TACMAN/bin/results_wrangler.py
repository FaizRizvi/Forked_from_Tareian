#!/usr/bin/env python

import pandas as pd
import glob
import os



class RESULTS:
    def __init__(self, dir_, name):
        self.column_types = {
            'BIN_ID': "category",
            'TF_NAME': "category"}

        self.directory = dir_
        self.name = name

    def get_matrix(self):
        os.chdir(self.directory)
        file_list = glob.glob("*_" +  + MODE + ".bed")
        df_list = []

        for i in file_list:
            df = pd.read_csv(i, sep="\t", 
                            header=None,
                            usecols=[6, 10],
                            dtype=self.column_types,
                            names=self.column_types.keys())
            df["SET_ID"] = i
            df_list.append(df)

        con = pd.concat(df_list)

        del df_list

        gb = con.groupby(["SET_ID", df.BIN_ID.astype(object), df.TF_NAME.astype(object)]).count()

        self.matrix = gb

    def get(self, x):
        print (self.x)

CHIP_RESULTS = RESULTS("/Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/output/intersection/CHIP/DHS/", "CHIP")

MOODS_RESULTS = RESULTS("/Users/caz3so/workspaces/thesis/bin/MOTIF_ANALYSIS/output/intersection/MOODS/DHS/", "MOODS")

dfcon = pd.concat([gp6, gpm], axis=1)
dfcon = dfcon.fillna(0)

dfcon = dfcon.apply(lambda x: [y if y <= 1 else 1 for y in x])

gpcon4 = dfcon.groupby(["P6", "CHIP"]).size()
gpcon4 = gpcon4.reset_index()
gpcon4["join"] = gpcon4["P6"].apply(str) + "_" + gpcon4["CHIP"].apply(str)

