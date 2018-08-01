cd /Users/caz3so/Dropbox/thesis/data/20180729_GM12878_DNASE_ENCODE/hg19

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt
import glob
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np

df = pd.read_csv("ENCODE_DNASE_counts.txt", sep = "\t")

sns.regplot(x=df["ENCFF000SLR_masked"], y=df["ENCFF000SLR_hg19"], fit_reg=False)
sns.jointplot(x=df["ENCFF000SLR_masked"], y=df["ENCFF000SLR_hg19"], kind='hex', marginal_kws=dict(bins=30, rug=True))
sns.jointplot(x=df["ENCFF000SLR_masked"], y=df["ENCFF000SLR_hg19"], kind='kde', color="grey", space=0)


sns.jointplot(x=df["ENCFF000SLR_masked"], y=df["ENCFF000SLR_hg19"], kind='scatter')

sns.jointplot(x=df["sepal_length"], y=df["sepal_width"], kind='scatter')
