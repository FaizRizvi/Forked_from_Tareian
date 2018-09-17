cd /Users/caz3so/Dropbox/thesis/data/20180729_GM12878_DNASE_ENCODE/hg19
cd /home/tacazares/Dropbox/thesis/data/20180729_GM12878_DNASE_ENCODE/hg19

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np

df = pd.read_csv("ENCODE_DNASE_counts.txt", header=0, sep = "\t")

df = df[~df.eq(0).any(1)]

df = df.dropna()
df['Variance'] = df.var(axis=1)

df = df[df.Variance != 0]

df = df.drop(columns="Variance")

fig = plt.figure(1)

z= linkage(df,method ="ward")

plt.subplot(421)
sns.regplot(x=df["ENCFF000SLR_masked"], y=df["ENCFF000SLR_hg19"], fit_reg=False)

plt.subplot(422)
sns.regplot(x=df["ENCFF000SLF_masked"], y=df["ENCFF000SLF_hg19"], fit_reg=False)

plt.subplot(423)
sns.regplot(x=df["ENCFF000SLG_masked"], y=df["ENCFF000SLG_hg19"], fit_reg=False)

plt.subplot(424)
sns.regplot(x=df["ENCFF000SLP_masked"], y=df["ENCFF000SLP_hg19"], fit_reg=False)

plt.subplot(425)
sns.regplot(x=df["ENCFF000SLL_masked"], y=df["ENCFF000SLL_hg19"], fit_reg=False)

plt.subplot(426)
sns.regplot(x=df["ENCFF001CUR_masked"], y=df["ENCFF001CUR_hg19"], fit_reg=False)

plt.subplot(427)
sns.regplot(x=df["ENCFF001CWQ_masked"], y=df["ENCFF001CWQ_hg19"], fit_reg=False)

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

sns.clustermap(df, z_score=0)
