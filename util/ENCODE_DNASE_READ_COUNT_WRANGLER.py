cd /Users/caz3so/Dropbox/thesis/data/20180729_GM12878_DNASE_ENCODE/hg19

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt
import glob

# The * is not a regex, it just means "match anything"
# This matches datafile-0.csv, datafile-1.csv, etc.
filenames_hg19 = glob.glob("*hg19.bed_read_length.bed")
filenames_masked =  glob.glob("*masked.bed_read_length.bed")

list_of_hg19_dfs = [pd.read_csv(filename, header=None, sep="\t", usecols=[6]) for filename in filenames_hg19]
list_of_masked_dfs = [pd.read_csv(filename, header=None, sep="\t", usecols=[6]) for filename in filenames_masked]

# zip loops through TWO THINGS AT ONCE
# so you're looking at dataframe #1 and filename #1
# then dataframe #2 and filename #2
# etc
# and assigning that filename as a new column in the dataframe

zippy1 = zip(list_of_hg19_dfs, filenames_hg19)
zippy2 = zip(list_of_masked_dfs, filenames_masked)

for dataframe, filename in zippy1:
    dataframe['filename'] = filename.replace("_hg19.bed_read_length.bed", "")

for dataframe, filename in zippy2:
    dataframe['filename'] = filename.replace("_masked.bed_read_length.bed", "")

result = pd.concat(list_of_hg19_dfs)
result["Genome"] = "hg19"
result_masked = pd.concat(list_of_masked_dfs)
result_masked["Genome"] = "masked"
results_total = pd.concat([result_masked, result])

flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]

results_total.groupby('filename')[6].sum()

y = sns.violinplot(x="filename" , y=6 , hue="Genome", hue_order=["hg19", "masked"], inner=None, data= results_total, split=True)
y.set(yscale="log")
y.set(ylabel="Length (BP)")
y.set(xlabel="Isoreplicate")
for item in y.get_xticklabels():
    item.set_rotation(60)
y.legend(bbox_to_anchor=(1, 1), loc=2)

result.to_csv("Salmon_combined_counts.csv", sep="\t")