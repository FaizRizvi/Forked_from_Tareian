#!/usr/bin/env

""""Principal Component Analysis"""
# Goal: I want to replicate Emily Miraldi's PCA code from Matlab in Python. 

"""Parameters"""
# imports the necessary libraries
import pandas as pd
import seaborn as sn
import numpy as np
import os
import errno
from scipy import stats
from matplotlib.mlab import PCA

output_dir_base_name = 'PCA'
input_gene_expression_file = '102717_DC_geo_expression_d100_p01_log2.txt'
dataset = 'myeloid'  # add a dataset name, this will be used in figure titles and output files

loadingOpt = 1     # Plot loadings?  1 --> yes, 0 --> no
totWeights = 75    # What number of top genes would you like to see? ~10k results in a blob.

gene_expression_cutoff = 1 # filter genes < minimum 
visualize_gene_expression_cutoff = 1 # visualize histogram of sample maximum 

fontSize = 12 # font size for figures    

dataScale = 'log2' # Set the data scaling options

if not os.path.exists(os.path.dirname(output_dir_base_name)): #if directory does not exist make it
    try:
        os.makedirs(os.path.dirname(output_dir_base_name))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

figinf = output_dir_base_name + dataset + '_' + dataScale  # Make a base name for output figures

plots = [1 1 1 1];
numweights = 10000;

"""Import Files and Obtain Parameters"""    
# load in the gene expression file
df = pd.read_csv(input_gene_expression_file, header=0, sep='\t', index_col=0)
df_mapped = pd.read_csv(input_gene_expression_file, header=0, sep='\t')

sampleNames = list(df.columns.values) # get sample names
totSamps = len(sampleNames) # get # of samples
genesc = df[df.columns[0]] # get gene names
ncounts = df[df.columns[1:]] # get gene expression values
[vars, obs] = np.shape(ncounts)

"""Analysis"""
# Apply gene expression cutoff
df['geneMaxes'] = df.max(axis=1)
df_cutoff = df.loc[df['geneMaxes'] >= gene_expression_cutoff]
df_cutoff = df.drop(columns='geneMaxes')

# find sum and calculate mean expression for each gene
df_p1 = df_cutoff + 1 # add a pseudocount and create new df

# Apply data normalization
if dataScale == 'log2':
    gene_mean = df_p1.mean(axis=1) # sum the matrix w/o genes in cutoff 
	df_mean_norm = df_p1.divide(gene_mean, axis=0)
	dataScaled = np.log2(df_p1)
    typename = 'log2(FC)'
elif dataScale == 'zscore':
	dataScaled = stats.zscore(df_cutoff)
    typename = 'z-score'
else:
    dataScaled = df_p1
    typename = 'Raw'

# make sure that totWeights isn't bigger than number of genes
#totGenes = dataScaled.count(axis=1)
#plotweights = 1:min(totWeights,size(dataScaled,1));

result = PCA(dataScaled) 
result.fracs

x_plot = []
y_plot = []
for item in result.Y:
 x_plot.append(item[0])
 y_plot.append(item[1])
 
plot_tuples = zip(x_plot[0:75], y_plot[0:75])
plot = pd.DataFrame(plot_tuples, columns = ("x","y"))

sn.regplot(x="x", y="y", data=plot, fit_reg=False, color='k')
#plt.ylim(10, -10)
#plt.xlim(-35, 35)
plt.axhline(y=0, c='silver')
plt.axvline(x=0, c='silver')
plt.marker

meta = pd.read_csv('/Users/caz3so/Dropbox/thesis/data/20180613_Myeloid_Figures/DC_GEOinfo_complete.txt', sep="\t", header=0)
df['Lineage'] = df_mapped['Name'].map(meta.set_index('Name')['Lineage'])

df = df.transpose()

# Separating out the features
x = df.loc[:, genesc].values
# Separating out the target
y = df.loc[:,['lineage']].values












