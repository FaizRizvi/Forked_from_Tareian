cd /Users/caz3so/Dropbox/thesis/data/201808012_ENCODE_DNASE_FINAL

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt
import matplotlib.pyplot as plt

# use the following code to generate the tables
#for i in *read_length.bed; do wc -l $i >> ENCODE_peaks.txt; done
# Create a dataframe with the four feature variables

df = pd.read_csv('ENCODE_PEAKS.txt', header = 0, sep='\t')
df
snsplot = sns.barplot(x="Sample", y="Peaks", hue="Library", data=df )
snsplot.legend(bbox_to_anchor=(1, 1), loc=2)

for item in snsplot.get_xticklabels():
    item.set_rotation(60)

pullsns = snsplot.get_figure()
pullsns.savefig("BARCHART_ENCODE.png", dpi=1000)

ax = sns.boxplot(x="Genome", y="Peaks", data=df)
ax = sns.swarmplot(x="Genome", y="Peaks", data=df, color=".25")
