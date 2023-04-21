#Train a classifier to predict the phenotypes from DESEQ genes dataset

#python3 multiclass_logistic_regression.py --path_input ../../data/HBC/counts_matrix_hvg.tab --path_label ../../data/HBC/cell_labels_oh.csv --path_test ../../data/HBC/counts_matrix_hvg.tab --path_save ../../data/HBC/hvg_MCLR_prob.csv

#Testing for multi-class -
#python multiclass_logistic_regression.py --path_input ../../data/HBC/counts_matrix_all.tab --path_label ../../data/HBC/cell_labels_oh.csv --path_test ../../data/CIBERSORTx/Expression_datasets/Fig2b-WholeBlood_RNAseq.hbc.all.tab --path_save ../../data/CIBERSORTx/Expression_datasets/MCLR_Fig2b-WholeBlood_RNAseq.hbc.all.csv
import numpy as np
import pandas as pd
import os
from numpy import savetxt
import random
import argparse

parser = argparse.ArgumentParser(description='Multi-class logistic regression')
parser.add_argument('--path_input', type=str, default='../../data/HBC/', help='path to training file')
parser.add_argument('--path_label', type=str, default='../../data/HBC/', help='path to label file')
parser.add_argument('--path_test', type=str, default='../../data/HBC/', help='path to test file')
parser.add_argument('--path_save', type=str, default='../../data/HBC/', help='path to results')

args = parser.parse_args()

with open(args.path_input,'r') as f:
    X = pd.read_csv(f,sep="\t")

with open(args.path_label,'r') as f:
    Y = pd.read_csv(f,sep=",",header=None)

with open(args.path_test,'r') as f:
    test_X = pd.read_csv(f,sep="\t")

label_all = np.where(Y==1)[1]
print(label_all)

data_X = X
data_X.drop(columns=data_X.columns[0], axis=1,  inplace=True)
#data_X.drop(data_X.index[0])

from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(data_X, label_all, test_size=0.20, random_state=0)
#np.savetxt("X_train.csv",x_train)
#np.savetxt("X_test.csv",x_test)
#np.savetxt("Y_train.csv",y_train)
#np.savetxt("Y_test.csv",y_test)

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
# all parameters not specified are set to their defaults
logisticRegr = LogisticRegression(multi_class='multinomial', solver='lbfgs',penalty='l2', C=1.0,max_iter=1000)
logisticRegr.fit(x_train, y_train)
print('data_X')
print(data_X)
predictions = logisticRegr.predict(x_test)

# Use score method to get accuracy of model
score = logisticRegr.score(x_test, y_test)
print(score)

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import metrics
cm = metrics.confusion_matrix(y_test, predictions)
print(cm)

plt.figure(figsize=(9,9))
sns.heatmap(cm, annot=True, fmt=".3f", linewidths=.5, square = True, cmap = 'Blues_r');
plt.ylabel('Actual label');
plt.xlabel('Predicted label');
all_sample_title = 'Accuracy Score: {0}'.format(score)
plt.title(all_sample_title, size = 15);

'''
print('Classes')
print(logisticRegr.classes_)
print('Coef')
print(logisticRegr.coef_)
np.savetxt("coef.csv",logisticRegr.coef_)
print('Intercept')
print(logisticRegr.intercept_)
'''

#Predictions for training set
test_X.drop(columns=test_X.columns[0], axis=1,  inplace=True)
predictions_proba = logisticRegr.predict_proba(test_X)
print('predict_probabilities')
print(predictions_proba)
np.savetxt(args.path_save,predictions_proba,delimiter=',')
