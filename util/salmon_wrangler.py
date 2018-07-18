cd /Users/caz3so/Dropbox/thesis/data/20180712_GEUVADIS_SALMON_QUANTS/quants

import glob
import seaborn as sns
import scipy.cluster.hierarchy as sch

# The * is not a regex, it just means "match anything"
# This matches datafile-0.csv, datafile-1.csv, etc.
filenames = glob.glob("*.sf")

list_of_dfs = [pd.read_csv(filename, header=0, sep="\t", usecols=[0,4], names=["genes", filename], skiprows=0, index_col=0) for filename in filenames]

# zip loops through TWO THINGS AT ONCE
# so you're looking at dataframe #1 and filename #1
# then dataframe #2 and filename #2
# etc
# and assigning that filename as a new column in the dataframe
#for dataframe, filename in zip(list_of_dfs, filenames):
#    dataframe['filename'] = filename

result = pd.concat(list_of_dfs, axis=1)
result['summation'] = result.sum(axis=1)
result = result[result.summation >= 10]
result = result.sort_values("summation")
result = result.drop("summation", axis=1)

top_100 = result.head(100)
bottom_100 = result.tail(100)
top_bottom = top_100.append(bottom_100)

heater = sns.heatmap(top_bottom, xticklabels=False, yticklabels=False)

result.to_csv("Salmon_combined_counts.csv", sep="\t")

Z = sch.linkage(result, method = 'ward')
den = sch.dendrogram(Z)
plt.title('Dendrogram for the clustering of the dataset on three different varieties of wheat kernels (Kama, Rosa and Canadian)')
plt.xlabel('Wheat kernels')
plt.ylabel('Euclidean distance in the space with dimensions AREA, PERIMETER and ASYMMETRY');
plt.show()

