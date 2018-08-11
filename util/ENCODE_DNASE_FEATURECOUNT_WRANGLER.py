cd /Users/caz3so/Dropbox/thesis/data/201808012_ENCODE_DNASE_FINAL/featureCounts

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob

filenames = glob.glob("*bam_counts.txt")

list_of_counts = [pd.read_csv(filename, header = 0, usecols=[6], sep="\t") for filename in filenames]
result = pd.concat(list_of_counts, axis=1)
result.to_csv('ENCODE_featureCounts.txt', sep='\t',index=False)

pairs = pd.read_csv("dnase_pairs.txt", sep= "\t")
zipper = zip(pairs.hg19, pairs.masked) 

result = pd.read_csv('ENCODE_featureCounts.txt', sep='\t')

fig = plt.figure()

for i, j in zipper:
    plt.figure()
    sns.kdeplot(result[i], shade=True)
    sns.kdeplot(result[j], shade=True)

result2 = result[(result != 0).all(1)]

result_heat = result

heater = sns.clustermap(result)
heater = sns.clustermap(result2, robust=True)

result2 = result.drop('Geneid', axis=1)
regro = sns.jointplot(x="ENCFF000SLL_hg19", y="ENCFF000SLL_masked", data= result, kind="kde")

g = sns.PairGrid(result)

sns.jointplot(x='ENCFF000SLL_hg19', y='ENCFF000SLL_hg19', data=result)

sns.pairplot(result)
y = sns.violinplot(x="filename" , y="Length" , hue="Genome", hue_order=["hg19", "masked"], inner=None, data= results_all, split=True)
y.set(yscale="log")
y.set(ylabel="Length (BP)")
y.set(xlabel="Sample")
for item in y.get_xticklabels():
    item.set_rotation(60)
y.legend(bbox_to_anchor=(1, 1), loc=2)


