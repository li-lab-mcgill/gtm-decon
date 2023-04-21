#Randomly shuffle the bulkRS data and split it into training and test sets (80:20)

#For TCGA-BRCA:
#python3 train_test_split_bulkRS_BRCA.py --path_input ~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.uniq.stage.tab --path_save ~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/stage/
#python3 train_test_split_bulkRS_BRCA.py --path_input ~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.hvg.TNBC_ER.tab --path_save ~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA_hvg_TNBC_ER/
#python3 train_test_split_bulkRS_BRCA.py --path_input ~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.hvg.type_2sets.tab --path_save ~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA_hvg_type_2sets/

import numpy as np
import pandas as pd
import os
from numpy import savetxt
import random
import argparse

parser = argparse.ArgumentParser(description='Train_test split datasets')
parser.add_argument('--path_input', type=str, default='../../../RA/other_datasets/HBC_NSR/', help='path to input files')
parser.add_argument('--path_save', type=str, default='../../../RA/other_datasets/HBC_NSR/train_valid_test/', help='path to save results')

args = parser.parse_args()

#~/projects/scDECON_CANCER/data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.stage.tab
with open(args.path_input,'r') as f:
        X = pd.read_csv(f,sep="\t")

orig_X = X

#Generate cell labels to one-hot encoding mapping
#Option 2: Label map
#label_map={'stage_i': 0, 'stage_ia': 1, 'stage_ib': 2, 'stage_ii': 3, 'stage_iia': 4, 'stage_iib': 5, 'stage_iii': 6, 'stage_iiia': 7, 'stage_iiib': 8, 'stage_iiic': 9, 'stage_iv': 10, 'stage_x': 11}
#label_map={'stage_i': 0, 'stage_ii': 1, 'stage_iii': 2}
#label_all =  np.array([label_map[i] for i in X['pathologic_stage']])

#label_map={'infiltrating_ductal_carcinoma': 0, 'infiltrating_lobular_carcinoma': 1, 'medullary_carcinoma': 2, 'metaplastic_carcinoma': 3, 'mixed_histology': 4, 'mucinous_carcinoma': 5, 'other': 6}
label_map={'infiltrating_ductal_carcinoma': 0, 'infiltrating_lobular_carcinoma': 1}
label_all =  np.array([label_map[i] for i in X['histological_type']])

#label_map={'TNBC': 0, 'ER+': 1}
#label_all =  np.array([label_map[i] for i in X['Subtype']])

#Save cell labels as one hot encoding
a = np.array(label_all)
b = np.zeros((a.size, a.max()+1))
b[np.arange(a.size),a] = 1
print('Saving cell labels...')
np.savetxt(args.path_save+'cell_labels_oh.csv',b,delimiter=",")

orig_X.drop('Unnamed:0',axis=1,inplace=True)
#orig_X.drop('pathologic_stage',axis=1,inplace=True)
orig_X.drop('histological_type',axis=1,inplace=True)
#orig_X.drop('Subtype',axis=1,inplace=True)
#orig_X.drop('PAM50',axis=1,inplace=True)
#orig_X.drop('ER_IHC',axis=1,inplace=True)
print('Saving counts matrix...')
orig_X.to_csv(args.path_save+'counts_matrix.tab',sep="\t")

X = orig_X
prior = pd.DataFrame(b)

print('shuffle')
print(len(orig_X))

index = list(range(len(orig_X)))
random.shuffle(index)

X_shuffled = orig_X.iloc[index]
prior_shuffled = prior.iloc[index]

'''
print(X_shuffled)
print(prior_shuffled)
print(Xn_shuffled)
print(Xln_shuffled)
print(Xv_shuffled)
print(Xvn_shuffled)
'''

#Split training:validation:test in 70:10:20 ratio
num_training = int(len(orig_X) * 0.8)

train_X = X_shuffled[:num_training]
train_X.drop(train_X.columns[0], axis = 1, inplace = True)
train_prior = prior_shuffled[:num_training]
train_prior = train_prior.to_numpy()

test_X = X_shuffled[num_training:]
test_X.drop(test_X.columns[0], axis = 1, inplace = True)
test_prior = prior_shuffled[num_training:]
test_prior = test_prior.to_numpy()

train_X.to_csv(args.path_save+'counts_matrix_train.tab',sep="\t")
savetxt(args.path_save+'cell_labels_oh_train.csv',train_prior,delimiter=",")
test_X.to_csv(args.path_save+'counts_matrix_test.tab',sep="\t")
savetxt(args.path_save+'cell_labels_oh_test.csv',test_prior,delimiter=",")

