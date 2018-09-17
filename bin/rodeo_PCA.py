
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import scale
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch

#import the csv and create a dataframe with 1st column as index, convert to NP array
data = pd.read_csv('102717_DC_geo_expression_d100_p01_log2.txt', sep='\t', header=0)
data = data.set_index('Name')
X=data.values

#calculate the d-dimensional mean vectors
mean_x = np.mean(X[0,:])
mean_y = np.mean(X[1,:])
mean_z = np.mean(X[2,:])
mean_vector = np.array([[mean_x],[mean_y],[mean_z]])

#create the covariance matrix for the data
cov_mat = np.cov([X[0,:],X[1,:],X[2,:]])
eig_val_cov, eig_vec_cov = np.linalg.eig(cov_mat)

#Sort eigenvectors and pairs
for ev in eig_vec_cov:
    np.testing.assert_array_almost_equal(1.0, np.linalg.norm(ev))

eig_pairs = [(np.abs(eig_val_cov[i]), eig_vec_cov[:,i]) for i in range(len(eig_val_cov))]
eig_pairs.sort(key=lambda x: x[0], reverse=True)

transformed = matrix_w.T.dot(X)

