cd /home/tacazares/Dropbox/thesis/data/201808012_ENCODE_DNASE_FINAL/featureCounts
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt
import glob

filenames = glob.glob("*_counts.txt")

list_of_counts = [pd.read_csv(filename, header = 0, usecols=[6], index_col = 0) for filename in filenames]
list_of_lengths = [pd.read_csv(filename, header = 0, usecols=[5], index_col = 0) for filename in filenames]

result = pd.concat(list_of_counts, axis=1)
result
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
