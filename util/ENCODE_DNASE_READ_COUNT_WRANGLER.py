cd /Users/caz3so/Dropbox/thesis/data/20180809_LCL_DNASE_MW_ENCODE/featurecounts//

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt
import glob

#Here we are creating a ilst of the filenames based on the source
filenames_hg19 = glob.glob("*_hg19.bed.bam_counts.txt_headless.txt")
filenames_masked =  glob.glob("*_masked.bed.bam_counts.txt_headless.txt")
filenames_MW = glob.glob("GM*.bed.bam_counts.txt_headless.txt")

#Here we are creating the dataframes based on the list of names
list_of_hg19_dfs = [pd.read_csv(filename, header=0, usecols=[0,5,6], sep="\t") for filename in filenames_hg19]
list_of_masked_dfs = [pd.read_csv(filename, header=0, usecols=[0,5,6], sep="\t") for filename in filenames_masked]
list_of_MW_dfs = [pd.read_csv(filename, header=0, usecols=[0,5,6], sep="\t") for filename in filenames_MW]

#We then want to create a zipped list of names and DFs and add a name for the source of the data
zippy1 = zip(list_of_hg19_dfs, filenames_hg19)
zippy2 = zip(list_of_masked_dfs, filenames_masked)
zippy3 = zip(list_of_MW_dfs, filenames_MW)

# Here we add the filename as an additional column to every dataframe and remove extra information
for dataframe, filename in zippy1:
    dataframe['filename'] = filename.replace("_hg19.bed.bam_counts.txt_headless.txt", "")

for dataframe, filename in zippy2:
    dataframe['filename'] = filename.replace("_masked.bed.bam_counts.txt_headless.txt", "")

for dataframe, filename in zippy3:
    dataframe['filename'] = filename.replace(".bed.bam_counts.txt_headless.txt", "")


# We then concatenate the files into one large frame for each source and add a column with the genome source
result = pd.concat(list_of_hg19_dfs)
result["Genome"] = "hg19"

result_masked = pd.concat(list_of_masked_dfs)
result_masked["Genome"] = "masked"

result_MW = pd.concat(list_of_MW_dfs)
result_MW["Genome"] = "MW"

results_total = pd.concat([result_masked, result, result_MW])

y = sns.violinplot(x="filename" , y="Length", inner=None, data= results_total)
y.set(yscale="log")
y.set(ylabel="Length (BP)")
y.set(xlabel="Sample")
for item in y.get_xticklabels():
    item.set_rotation(60)
y.legend(bbox_to_anchor=(1, 1), loc=2)

