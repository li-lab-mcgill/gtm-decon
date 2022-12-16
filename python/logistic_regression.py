#Train a classifier to predict the phenotypes from DESEQ genes dataset

#python3 logistic_regression.py --path_input ~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA_TNBC_ER_DESEQ2.tab
#python3 logistic_regression.py --path_input ~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.TNBC_ER.tab
#python3 logistic_regression.py --path_input ~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.TNBC_ER.sparse.75p.tab

#python3 logistic_regression.py --path_input /home/mcb/users/slaksh1/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA_type_lobular_vs_ductal_DESEQ2.tab
#python3 logistic_regression.py --path_input /home/mcb/users/slaksh1/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.uniq.type.2sets.tab
import numpy as np
import pandas as pd
import os
from numpy import savetxt
import random
import argparse

parser = argparse.ArgumentParser(description='Logistic regression')
parser.add_argument('--path_input', type=str, default='../../../RA/other_datasets/HBC_NSR/', help='path to input files')

args = parser.parse_args()

#~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.stage.tab
with open(args.path_input,'r') as f:
    X = pd.read_csv(f,sep="\t")


#TNBC-ER+
#subtype = np.array(X['Subtype'])
#label_map={'TNBC': 0, 'ER+': 1}
#label_all =  np.array([label_map[i] for i in X['Subtype']])

#data_X = X
#data_X.drop('Unnamed:0',axis=1,inplace=True)
#data_X.drop('Subtype',axis=1,inplace=True)

#data_X.drop('PAM50',axis=1,inplace=True)
#data_X.drop('ER_IHC',axis=1,inplace=True)

#ductal carcinoma
subtype = np.array(X['histological_type'])
data_X = X
data_X.drop('Unnamed:0',axis=1,inplace=True)
data_X.drop('histological_type',axis=1,inplace=True)

from sklearn.model_selection import train_test_split
x_train, x_test, y_train, y_test = train_test_split(data_X, subtype, test_size=0.20, random_state=0)

from sklearn.linear_model import LogisticRegression
# all parameters not specified are set to their defaults
logisticRegr = LogisticRegression(max_iter=1000)
logisticRegr.fit(x_train, y_train)
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
