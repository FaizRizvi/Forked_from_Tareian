#!/usr/bin/env python

import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import seaborn as sns
import pandas as pd

"""Preferences"""
colors = sns.color_palette()
csv_delimiter = '\t'
index_col_number = 0
file_name = 'jake_data.txt'
group_number = [0, 1, 2, 3, 4, 5, 6, 7]

"""Database creation and manipulation"""
#create dataframe
data = pd.read_csv(file_name,
                sep = csv_delimiter,
                header = 0,
                index_col=index_col_number)

#you must transpose the data in order to have the PCA work correctly. Columns= genes and rows=samples
transposed_data = data.transpose()

"""Set up name/number/colors"""
#you can change the target names by adding in an array of values. You might need to add a column to DF to track the names with the samples.
index_names = list(transposed_data.index.values)
target_names = index_names
color_names = zip(colors, group_number, target_names)

"""PCA"""
pca = PCA()

pca.fit(transposed_data)

transformed_data = pca.transform(transposed_data)

pc1f=float((pca.explained_variance_ratio_[0]))
pc2f=float((pca.explained_variance_ratio_[1]))

pc1 = str(100*round(pc1f, 2))
pc2 = str(100*round(pc2f, 2))

"""Plotting"""
#get explained variance and convert to float
pc1f=float((pca.explained_variance_ratio_[0]))
pc2f=float((pca.explained_variance_ratio_[1]))

pc1 = str(100*round(pc1f, 2))
pc2 = str(100*round(pc2f, 2))

with plt.style.context('seaborn-whitegrid'):
    plt.figure()
    for c, i, labels in color_names:
        plt.scatter(transformed_data[i:, 0],
                    transformed_data[i:,1],
                    label=labels)
    #add labels
    plt.xlabel('PC1: ' + pc1 + '%')
    plt.ylabel('PC2: ' + pc2 + '%')
    #add title
    plt.title('PCA of Jakes Data')
    #add legend
    plt.legend()
    #add verticle line
    plt.axvline(0, color='k')
    plt.axhline(0, color='k')
    #remove tick labels
    plt.axes().get_xaxis().set_ticks([])
    plt.axes().get_yaxis().set_ticks([])

plt.show()
