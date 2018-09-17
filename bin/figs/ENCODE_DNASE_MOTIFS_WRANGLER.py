cd /Users/caz3so/Dropbox/thesis/data/20180812_LCL_Dnase_hg19/

import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("GM12878_hg19_UNION_MOTIFS_SUBS.bed", header = None, sep="\t", names=["Chr", "Start", "Stop", "Subcompartment", "Motifs"])

df["Subcompartment"].fillna("Unknown", inplace=True)

df = df.sort_values("Subcompartment")

df

df["Length"] = df["Stop"] - df["Start"]

df['Motifs_count'] = df.Motifs.map(lambda x: [i.strip() for i in x.split(",")])

df['Motifs_count'] = df['Motifs'].str.len()

df["one"] = 1

equiv = {"A1":400600000, "A2":581400000, "B1":349400000, "B2":435700000, "B3":855500000, "B4":11000000, "Unknown":248700000}

group_counts = df.groupby(["Subcompartment"])
group_counts_sum["Sub_length"] = group_counts_sum["Subcompartment"].map(equiv)
group_counts_sum["Motifs/BP"] = group_counts_sum["Motifs_count"] / group_counts_sum["Sub_length"]
group_counts_sum["Motifs/Peaks/BP"] = (group_counts_sum["Motifs_count"] / group_counts_sum["one"] ) / group_counts_sum["Sub_length"]

#Barplot
t = sns.barplot(x="Subcompartment", y = "Motifs/Peaks/BP", data=group_counts_sum, palette=pink)
t.set(xlabel="Sub-compartment")
t.set(ylabel="TF Motifs/Peaks/BP")
t.set(title="TF Motifs/Peaks/BP by Sub-compartment")

#Countplot
y = sns.countplot(x="Subcompartment", data=union_frame, order=["A1", "A2", "B1", "B2", "B3", "B4"], palette=flatui)
for p in y.patches:
        y.annotate('{:}'.format(p.get_height()), (p.get_x(), p.get_height()))

y = sns.countplot(x="Subcompartment", data=df, palette=green)
for p in y.patches:
        y.annotate('{:}'.format(p.get_height()), (p.get_x(), p.get_height()))
y.set(ylabel="# of Peaks")
y.set(title="Peaks per Subcompartment in Union")

#Colors
flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
colors = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
colorful = sns.xkcd_palette(colors)
grade = sns.light_palette("green")
diverging = sns.diverging_palette(10, 220, sep=80, n=7)
green = sns.cubehelix_palette(8, start=.5, rot=-.75)
cube = sns.color_palette("cubehelix", 8)
pink = sns.cubehelix_palette(8)
teal = sns.color_palette("GnBu_d")
husl = sns.color_palette("husl", 8)

group_counts_sum.to_csv("groupby_sums.txt")
