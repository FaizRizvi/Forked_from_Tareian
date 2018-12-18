cd /Users/caz3so/Dropbox/Data/Roskin/data/IgH

# Load the library with the iris dataset
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score

# Create a random forest Classifier. By convention, clf means 'Classifier'
clf = RandomForestClassifier(n_estimators=200, n_jobs=-1)

#column types
column_types = {'specimen_tissue':        object,
                'cdr3_length':            np.int32,
                'cdr3_hydrophobicity':    np.float32,
                'cdr3_isoelectric_point': np.float32,
                'cdr3_any_nglycosylation':np.float32,
                'cdr3_instability':       np.float32,
                'cdr3_aromaticity':       np.float32,
                'cdr3_mut_freq':       np.float32}

# Create a dataframe with the four feature variables
df = pd.read_csv('BCR_randomforest.csv', header=0, sep=',', usecols=column_types.keys(), dtype=column_types)

features = ['cdr3_length', 'cdr3_hydrophobicity', 'cdr3_isoelectric_point', 'cdr3_any_nglycosylation', 'cdr3_instability', 'cdr3_aromaticity', 'cdr3_mut_freq']

#Fill Na
df = df.dropna()
df = df[df.specimen_tissue != '??']

# Create a new column that for each row, generates a random number between 0 and 1, and
# be used as the training data and some as the test data.
df['is_train'] = np.random.uniform(0, 1, len(df)) <= .75

# Create two new dataframes, one with the training rows, one with the test rows
train, test = df[df['is_train']], df[~df['is_train']==False]

# train['species'] contains the actual species names. Before we can use it convert to number
y, uniques = pd.factorize(train['specimen_tissue'])

# Train the Classifier to take the training features and learn how they relate
clf.fit(train[features], y)

#apply prediction to test set
preds = clf.predict(test[features])
uniques = uniques[preds]
# Create confusion matrix
pd.crosstab(test['specimen_tissue'], preds, rownames=['Actual'], colnames=['Predicted'])

# View a list of the features and their importance scores
list(zip(train[features], clf.feature_importances_))

# View The Accuracy Of Our Full Feature (4 Features) Model
print "Train Accuracy :: ", accuracy_score(test['specimen_tissue'], preds)
