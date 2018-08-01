import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt

# use the following code to generate the tables
#for i in *read_length.bed; do wc -l $i >> ENCODE_peaks.txt; done
# Create a dataframe with the four feature variables

df = pd.read_csv('ENCODE_peaks.txt', header = 0,sep='\t')

ax = sns.barplot(x="Isoreplicate", y="Peaks", hue="Genome", data=df )
ax.legend(bbox_to_anchor=(1, 1), loc=2)

for item in ax.get_xticklabels():
    item.set_rotation(60)
 
 
ax = sns.boxplot(x="Genome", y="Peaks", data=df)
ax = sns.swarmplot(x="Genome", y="Peaks", data=df, color=".25")

