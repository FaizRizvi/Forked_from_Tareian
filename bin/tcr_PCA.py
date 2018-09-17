cd /Users/caz3so/Dropbox/Data/Roskin/data/TCRB

#!/usr/bin/env python

import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import seaborn as sns
import pandas as pd

"""Database creation and manipulation"""
#create dataframe
df = pd.read_csv('TCR_clone_v_tissue_size_population_alleles.csv',
                sep=',',
                header=0,
                usecols= ('v_segment', 'specimen_tissue', 'v_segment.size.size'))
             
#Create pivot table
df2 = df.pivot(index='v_segment',
                columns='specimen_tissue',
                values='v_segment.size.size')
                
#Fill na
df2 = df2.fillna(value=0)

#transpose
df2 = df2.transpose()
#drop ?? tissues
df2 = df2.drop('??')

"""Set up name/number/colors"""
colors = sns.color_palette()
target_names = list(df2.index.values)
number = [0, 1, 2, 3, 4, 5, 6, 7]
color_names = zip(colors, number, target_names)

"""PCA"""
#standardize data
X_std = StandardScaler().fit_transform(df2)

#declare number of components
sklearn_pca = PCA(n_components=1)

#PCAand transform data
Y_sklearn = sklearn_pca.fit_transform(X_std)

"""Plotting"""
#get explained variance and convert to float
pc1f=float((sklearn_pca.explained_variance_[0]))
pc2f=float((sklearn_pca.explained_variance_[1]))

#take float and roud to two decimels and then convert to string
pc1 = str(round(pc1f, 2))
pc2 = str(round(pc2f, 2))

#plot each tissue as a different color using seaborn
with plt.style.context('seaborn-whitegrid'):
    plt.figure()
    for c, i, labels in color_names:
        plt.scatter(Y_sklearn[i:, 0],
                    Y_sklearn[i:,1],
                    label=labels)
                    
    #add labels
    plt.xlabel('PC1: ' + pc1 + '%')
    plt.ylabel('PC2: ' + pc2 + '%')
    
    #add title
    plt.title('PCA of TCR V-segment usage individual')
    
    #add legend
    plt.legend()
    
    #add verticle line
    plt.axvline(0, color='k')
    plt.axhline(0, color='k')
    
    #remove tick labels
    plt.axes().get_xaxis().set_ticks([])
    plt.axes().get_yaxis().set_ticks([])
    
plt.show()

pca = PCA().fit(df2)
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance');

