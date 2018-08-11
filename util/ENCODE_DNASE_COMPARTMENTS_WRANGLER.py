cd /Users/caz3so/Dropbox/thesis/data/20180809_LCL_DNASE_MW_ENCODE/BED/format_4col/

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob

#Here we are creating a ilst of the filenames based on the source
filenames_hg19 = glob.glob("*_hg19.bed_subs.bed")
filenames_masked =  glob.glob("*_masked.bed_subs.bed")
filenames_MW = glob.glob("GM*.bed_subs.bed")

#Here we are creating the dataframes based on the list of names
list_of_hg19_dfs = [pd.read_csv(filename, header=None, sep="\t") for filename in filenames_hg19]
list_of_masked_dfs = [pd.read_csv(filename, header=None, sep="\t") for filename in filenames_masked]
list_of_MW_dfs = [pd.read_csv(filename, header=None, sep="\t") for filename in filenames_MW]

#We then want to create a zipped list of names and DFs and add a name for the source of the data
zippy1 = zip(list_of_hg19_dfs, filenames_hg19)
zippy2 = zip(list_of_masked_dfs, filenames_masked)
zippy3 = zip(list_of_MW_dfs, filenames_MW)

# Here we add the filename as an additional column to every dataframe and remove extra information
for dataframe, filename in zippy1:
    dataframe['filename'] = filename.replace("_hg19.bed_subs.bed", "")

for dataframe, filename in zippy2:
    dataframe['filename'] = filename.replace("_masked.bed_subs.bed", "")

for dataframe, filename in zippy3:
    dataframe['filename'] = filename.replace(".bed_subs.bed", "")

# We then concatenate the files into one large frame for each source and add a column with the genome source
result = pd.concat(list_of_hg19_dfs)
result["Genome"] = "hg19"

result_masked = pd.concat(list_of_masked_dfs)
result_masked["Genome"] = "masked"

result_MW = pd.concat(list_of_MW_dfs)
result_MW["Genome"] = "MW"

#Here we combine the indiviual files into genome level files 
results_total = pd.concat([result_masked, result, result_MW])

#Change the column names/add them
results_total.columns = ["Chr", "Start", "Stop", "Subcompartment", "some", "random", "start", "stop", "color", "filename", "Genome"]

#get the lengths of the peaks and add them as a new row
results_total["Length"] = results_total["Stop"] - results_total["Start"]

#change the Nan values in the sub column to unknown
results_total["Subcompartment"].fillna("Unknown", inplace=True)

#create a dict of sub:length keys- made from summing and hard coding
equiv = {"A1":400600000, "A2":581400000, "B1":349400000, "B2":435700000, "B3":855500000, "B4":11000000, "Unknown":248700000}

#Get the grouby object
group_counts = results_total.groupby(["Genome","filename","Subcompartment"])

#create an object of sizes of subcompartments next to the number of peaks in each category
group_counts_size = group_counts.size()
group_counts_size_frame = group_counts_size.to_frame()
group_counts_size_frame = group_counts_size.reset_index()

#map the normalized lengths from a reference
group_counts_size_frame["Sub_length"] = group_counts_size_frame["Subcompartment"].map(equiv)
group_counts_size_frame["Peaks/BP"] = group_counts_size_frame[0] / group_counts_size_frame["Sub_length"]
group_counts_size_frame.columns = ["Genome", "filename", "Subcompartment", "Peaks", "Sub_length", "Peaks/BP"]

#sums
group_counts_sum = group_counts.sum()
group_counts_sum_frame = group_counts_sum.to_frame()
group_counts_sum = group_counts_sum.reset_index()

#Colors
flatui = ["#ef8a62", "#ef8a62", "#67a9cf", "#67a9cf", "#67a9cf", "#67a9cf", "#ffffbf"]
colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
colorful = sns.xkcd_palette(colors)
grade = sns.light_palette("green")
diverging = sns.diverging_palette(10, 220, sep=80, n=7)

#Boxplot
t = sns.boxplot(x="Subcompartment", y = "Peaks/BP", hue="Genome", hue_order=["hg19"], data=group_counts_size_frame, palette=flatui)
plt.yscale("log")

#Countplot
y = sns.countplot(x="Subcompartment", data=union_frame, order=["A1", "A2", "B1", "B2", "B3", "B4"], palette=flatui)
for p in y.patches:
        y.annotate('{:}'.format(p.get_height()), (p.get_x(), p.get_height()))

group_counts_sum
