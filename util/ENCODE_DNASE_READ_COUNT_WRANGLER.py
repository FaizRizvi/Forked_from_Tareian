cd /Users/caz3so/Dropbox/thesis/data/20180729_GM12878_DNASE_ENCODE/hg19

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt
import glob

filenames_hg19 = glob.glob("*hg19.bed_read_length.bed")
filenames_masked =  glob.glob("*masked.bed_read_length.bed")

list_of_hg19_dfs = [pd.read_csv(filename, header=None, sep="\t", usecols=[6]) for filename in filenames_hg19]
list_of_masked_dfs = [pd.read_csv(filename, header=None, sep="\t", usecols=[6]) for filename in filenames_masked]

df = pd.read_csv('ENCODE_peak_union_ID_lengths.txt', header=0, sep="\t", usecols=[5])

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

df["Genome"] = "hg19"
df["filename"] = "Union"

results_all = pd.concat([results_total, df])

y = sns.violinplot(x="filename" , y="Length" , hue="Genome", hue_order=["hg19", "masked"], inner=None, data= results_all, split=True)
y.set(yscale="log")
y.set(ylabel="Length (BP)")
y.set(xlabel="Sample")
for item in y.get_xticklabels():
    item.set_rotation(60)
y.legend(bbox_to_anchor=(1, 1), loc=2)
