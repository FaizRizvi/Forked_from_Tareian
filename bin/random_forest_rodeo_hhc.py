cd /Users/caz3so/Dropbox/Data/Roskin/data/IgH/
cd /Users/caz3so/Dropbox/Data/Roskin/data/HHC/individual/
cd /Users/caz3so/Dropbox/Data/Roskin/data/HHC/individual2/

# Load the library with the iris dataset
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn import tree
import graphviz
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

features = ['cdr3_hydrophobicity', 'cdr3_isoelectric_point', 'cdr3_instability']
coltouse = ['specimen_tissue', 'cdr3_hydrophobicity', 'cdr3_isoelectric_point', 'cdr3_instability']

#column types
column_types = {'specimen_tissue':        object,
                'cdr3_length':            np.int32,
                'cdr3_hydrophobicity':    np.float32,
                'cdr3_isoelectric_point': np.float32,
                'cdr3_any_nglycosylation':np.float32,
                'cdr3_instability':       np.float32,
                'cdr3_aromaticity':       np.float32,
                'cdr3_mut_freq':       np.float32}

df = pd.read_csv('BCR_randomforest.csv', header=0, sep=',', usecols=coltouse, dtype=column_types)
df2 = pd.read_csv('PBMC_HHC.csv', header=0, sep=',', usecols=coltouse, dtype=column_types)

#Fill Na
df = df.dropna()
df = df[df.specimen_tissue != '??']

df2 = df2.dropna()
df2 = df2[df2.specimen_tissue != '??']

df2['specimen_tissue'] = 'PBMC-HHC'

df3 = pd.concat([df, df2])
df4 = pd.read_csv('PBMC_HHC2.csv', header=0, sep=',', usecols=coltouse, dtype=column_types)
df4 = df4.dropna()
df4 = df4[df4.specimen_tissue != '??']
df4['specimen_tissue'] = 'PBMC-HHC'

# are three species, which have been coded as 0, 1, or 2.
y, names = pd.factorize(df3['specimen_tissue'])

# Create a random forest Classifier. By convention, clf means 'Classifier'
clf = RandomForestClassifier(n_estimators=100, n_jobs=-1)

# Train the Classifier to take the training features and learn how they relate
# to the training y (the species)
clf.fit(df3[features], y)

#apply prediction to test set
preds = clf.predict(df4[features])
trainpreds = clf.predict(df3[features])

# Create confusion matrix
cross = pd.crosstab(df4['specimen_tissue'], names[preds], rownames=['Actual'], colnames=['Predicted'])

# View a list of the features and their importance scores
zipped = list(zip(df4[features], clf.feature_importances_))

# View The Accuracy Of Our Full Feature (4 Features) Model
print "Test Accuracy :: ", accuracy_score(df4['specimen_tissue'], names[preds])
print "Train Accuracy :: ", accuracy_score(df3['specimen_tissue'], names[trainpreds])

clf = tree.DecisionTreeClassifier()
clf = clf.fit(df[features], y)
clf.predict([[2., 2.]])
dot_data = tree.export_graphviz(clf, out_file=None)
graph = graphviz.Source(dot_data)
graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
graph = pydotplus.graph_from_dot_data(dot_data)
graph.write_png('tree.png')

plt.figure()
sns.boxplot(x="specimen_tissue", y="cdr3_instability", data=df3);
plt.show()
