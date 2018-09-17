cd /Users/caz3so/Dropbox/Data/Roskin/data/TCRB/

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

#column types
column_types = {'specimen_tissue':        object,
                'cdr3_length':            np.int32,
                'cdr3_hydrophobicity':    np.float32,
                'cdr3_isoelectric_point': np.float32,
                'cdr3_any_nglycosylation':np.float32,
                'cdr3_instability':       np.float32,
                'cdr3_aromaticity':       np.float32,
                'cdr3_mut_freq':       np.float32}

df = pd.read_csv('TCR_randomforest.csv', header=0, sep=',', usecols=column_types.keys(), dtype=column_types)

features = ['cdr3_hydrophobicity', 'cdr3_isoelectric_point', 'cdr3_instability']

#Fill Na
df = df.dropna()
df = df[df.specimen_tissue != '??']

# Create a new column that for each row, generates a random number between 0 and 1, and
# if that value is less than or equal to .75, then sets the value of that cell as True
# and false otherwise. This is a quick and dirty way of randomly assigning some rows to
# be used as the training data and some as the test data.
df['is_train'] = np.random.uniform(0, 1, len(df)) <= .75

# Create two new dataframes, one with the training rows, one with the test rows
train, test = df[df['is_train']], df[~df['is_train']]

# train['species'] contains the actual species names. Before we can use it,
# we need to convert each species name into a digit. So, in this case there
# are three species, which have been coded as 0, 1, or 2.
y, names = pd.factorize(train['specimen_tissue'])

# Create a random forest Classifier. By convention, clf means 'Classifier'
clf = RandomForestClassifier(n_estimators=100, n_jobs=-1)

# Train the Classifier to take the training features and learn how they relate
# to the training y (the species)
clf =clf.fit(train[features], y)

#apply prediction to test set
preds = clf.predict(test[features])
trainpreds = clf.predict(train[features])

# Create confusion matrix
pd.crosstab(test['specimen_tissue'], names[preds], rownames=['Actual'], colnames=['Predicted'])

# View a list of the features and their importance scores
list(zip(train[features], clf.feature_importances_))

# View The Accuracy Of Our Full Feature (4 Features) Model
print "Test Accuracy :: ", accuracy_score(test['specimen_tissue'], names[preds])
print "Train Accuracy :: ", accuracy_score(train['specimen_tissue'], names[trainpreds])

clf = tree.DecisionTreeClassifier()
clf = clf.fit(df[features], y)
clf.predict([[2., 2.]])
dot_data = tree.export_graphviz(clf, out_file=None)
graph = graphviz.Source(dot_data)
graph = pydotplus.graph_from_dot_data(dot_data.getvalue())
graph = pydotplus.graph_from_dot_data(dot_data)
graph.write_png('tree.png')


with open("fruit_classifier.txt", "w") as f:
    f = tree.export_graphviz(clf, out_file=f)
