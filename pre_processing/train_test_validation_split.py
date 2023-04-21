#Randomly shuffle the single cell data and split it into training and test sets

#For cases where diff number of samples are used as training etc.
#python3 train_test_validation_split.py --path_input ../../data/PI_Segerstolpe_scRS/all_patients/4_ --path_save ../../data/PI_Segerstolpe_scRS/all_patients/train_valid_test/4_
#python3 train_test_validation_split.py --path_input ../../data/PI_Segerstolpe_scRS/all_patients/6_ --path_save ../../data/PI_Segerstolpe_scRS/all_patients/train_valid_test/6_

import numpy as np
import pandas as pd
import os
from numpy import savetxt
import random
import argparse

parser = argparse.ArgumentParser(description='Train_test_validation split datasets')
parser.add_argument('--path_input', type=str, default='../../../RA/other_datasets/HBC_NSR/', help='path to input files')
parser.add_argument('--path_save', type=str, default='../../../RA/other_datasets/HBC_NSR/train_valid_test/', help='path to save results')

args = parser.parse_args()



# Get prior probabilities
with open(args.path_input+'cell_labels_oh.csv','r') as f:
    prior = pd.read_csv(f,sep=",",header=None)
#with open(args.path_input+'cell_labels_oh_7celltypes.csv','r') as f:
#    prior = pd.read_csv(f,sep=",",header=None)
#with open(args.path_input+'cell_labels_oh_ordered.csv','r') as f:
#    prior = pd.read_csv(f,sep=",",header=None)

with open(args.path_input+'counts_matrix_all.tab','r') as f:
    XA = pd.read_csv(f,sep="\t")

with open(args.path_input+'counts_matrix_all_normr.tab','r') as f:
    XAn = pd.read_csv(f,sep="\t")

with open(args.path_input+'counts_matrix_all_log1p_normr.tab','r') as f:
    XAln = pd.read_csv(f,sep="\t")

with open(args.path_input+'counts_matrix_pp.tab','r') as f:
    X = pd.read_csv(f,sep="\t")

with open(args.path_input+'counts_matrix_pp_normr.tab','r') as f:
    Xn = pd.read_csv(f,sep="\t")

with open(args.path_input+'counts_matrix_pp_log1p_normr.tab','r') as f:
    Xln = pd.read_csv(f,sep="\t")

with open(args.path_input+'counts_matrix_hvg.tab','r') as f:
    Xv = pd.read_csv(f,sep="\t")

with open(args.path_input+'counts_matrix_hvg_log1p_normr.tab','r') as f:
    Xvn = pd.read_csv(f,sep="\t")

print('shuffle')
print(len(X))

index = list(range(len(X)))
random.seed(1000)
random.shuffle(index)
#with open(args.path_input+'cell_labels_oh_ordered.csv','r') as f:
#    prior = pd.read_csv(f,sep=",",header=None)

X_shuffled = X.iloc[index]
prior_shuffled = prior.iloc[index]
Xn_shuffled = Xn.iloc[index]
Xln_shuffled = Xln.iloc[index]
Xv_shuffled = Xv.iloc[index]
Xvn_shuffled = Xvn.iloc[index]
XA_shuffled = XA.iloc[index]
XAn_shuffled = XAn.iloc[index]
XAln_shuffled = XAln.iloc[index]

'''
print(X_shuffled)
print(prior_shuffled)
print(Xn_shuffled)
print(Xln_shuffled)
print(Xv_shuffled)
print(Xvn_shuffled)
'''

#Split training:validation:test in 70:10:20 ratio
num_training = int(len(X) * 0.7)
num_validation = int(len(X) * 0.8)

train_X = X_shuffled[:num_training]
train_X.drop(train_X.columns[0], axis = 1, inplace = True)
train_prior = prior_shuffled[:num_training]
train_prior = train_prior.to_numpy()
train_Xn = Xn_shuffled[:num_training]
train_Xn.drop(train_Xn.columns[0], axis = 1, inplace = True)
train_Xln = Xln_shuffled[:num_training]
train_Xln.drop(train_Xln.columns[0], axis = 1, inplace = True)
train_Xv = Xv_shuffled[:num_training]
train_Xv.drop(train_Xv.columns[0], axis = 1, inplace = True)
train_Xvn = Xvn_shuffled[:num_training]
train_Xvn.drop(train_Xvn.columns[0], axis = 1, inplace = True)
train_XA = XA_shuffled[:num_training]
train_XA.drop(train_XA.columns[0], axis = 1, inplace = True)
train_XAn = XAn_shuffled[:num_training]
train_XAn.drop(train_XAn.columns[0], axis = 1, inplace = True)
train_XAln = XAln_shuffled[:num_training]
train_XAln.drop(train_XAln.columns[0], axis = 1, inplace = True)

validation_X = X_shuffled[num_training:num_validation]
validation_X.drop(validation_X.columns[0], axis = 1, inplace = True)
validation_prior = prior_shuffled[num_training:num_validation]
validation_prior = validation_prior.to_numpy()
validation_Xn = Xn_shuffled[num_training:num_validation]
validation_Xn.drop(validation_Xn.columns[0], axis = 1, inplace = True)
validation_Xln = Xln_shuffled[num_training:num_validation]
validation_Xln.drop(validation_Xln.columns[0], axis = 1, inplace = True)
validation_Xv = Xv_shuffled[num_training:num_validation]
validation_Xv.drop(validation_Xv.columns[0], axis = 1, inplace = True)
validation_Xvn = Xvn_shuffled[num_training:num_validation]
validation_Xvn.drop(validation_Xvn.columns[0], axis = 1, inplace = True)
validation_XA = XA_shuffled[num_training:num_validation]
validation_XA.drop(validation_XA.columns[0], axis = 1, inplace = True)
validation_XAn = XAn_shuffled[num_training:num_validation]
validation_XAn.drop(validation_XAn.columns[0], axis = 1, inplace = True)
validation_XAln = XAln_shuffled[num_training:num_validation]
validation_XAln.drop(validation_XAln.columns[0], axis = 1, inplace = True)

test_X = X_shuffled[num_validation:]
test_X.drop(test_X.columns[0], axis = 1, inplace = True)
test_prior = prior_shuffled[num_validation:]
test_prior = test_prior.to_numpy()
test_Xn = Xn_shuffled[num_validation:]
test_Xn.drop(test_Xn.columns[0], axis = 1, inplace = True)
test_Xln = Xln_shuffled[num_validation:]
test_Xln.drop(test_Xln.columns[0], axis = 1, inplace = True)
test_Xv = Xv_shuffled[num_validation:]
test_Xv.drop(test_Xv.columns[0], axis = 1, inplace = True)
test_Xvn = Xvn_shuffled[num_validation:]
test_Xvn.drop(test_Xvn.columns[0], axis = 1, inplace = True)
test_XA = XA_shuffled[num_validation:]
test_XA.drop(test_XA.columns[0], axis = 1, inplace = True)
test_XAn = XAn_shuffled[num_validation:]
test_XAn.drop(test_XAn.columns[0], axis = 1, inplace = True)
test_XAln = XAln_shuffled[num_validation:]
test_XAln.drop(test_XAln.columns[0], axis = 1, inplace = True)

train_X.to_csv(args.path_save+'counts_matrix_pp_train.tab',sep="\t")
savetxt(args.path_save+'cell_labels_oh_train.csv',train_prior,delimiter=",")
validation_X.to_csv(args.path_save+'counts_matrix_pp_validation.tab',sep="\t")
savetxt(args.path_save+'cell_labels_oh_validation.csv',validation_prior,delimiter=",")
test_X.to_csv(args.path_save+'counts_matrix_pp_test.tab',sep="\t")
savetxt(args.path_save+'cell_labels_oh_test.csv',test_prior,delimiter=",")

train_Xn.to_csv(args.path_save+'counts_matrix_pp_normr_train.tab',sep="\t")
validation_Xn.to_csv(args.path_save+'counts_matrix_pp_normr_validation.tab',sep="\t")
test_Xn.to_csv(args.path_save+'counts_matrix_pp_normr_test.tab',sep="\t")

train_Xln.to_csv(args.path_save+'counts_matrix_pp_log1p_normr_train.tab',sep="\t")
validation_Xln.to_csv(args.path_save+'counts_matrix_pp_log1p_normr_validation.tab',sep="\t")
test_Xln.to_csv(args.path_save+'counts_matrix_pp_log1p_normr_test.tab',sep="\t")

train_Xv.to_csv(args.path_save+'counts_matrix_hvg_train.tab',sep="\t")
validation_Xv.to_csv(args.path_save+'counts_matrix_hvg_validation.tab',sep="\t")
test_Xv.to_csv(args.path_save+'counts_matrix_hvg_test.tab',sep="\t")

train_Xvn.to_csv(args.path_save+'counts_matrix_hvg_log1p_normr_train.tab',sep="\t")
validation_Xvn.to_csv(args.path_save+'counts_matrix_hvg_log1p_normr_validation.tab',sep="\t")
test_Xvn.to_csv(args.path_save+'counts_matrix_hvg_log1p_normr_test.tab',sep="\t")

train_XA.to_csv(args.path_save+'counts_matrix_all_train.tab',sep="\t")
validation_XA.to_csv(args.path_save+'counts_matrix_all_validation.tab',sep="\t")
test_XA.to_csv(args.path_save+'counts_matrix_all_test.tab',sep="\t")

train_XAn.to_csv(args.path_save+'counts_matrix_all_normr_train.tab',sep="\t")
validation_XAn.to_csv(args.path_save+'counts_matrix_all_normr_validation.tab',sep="\t")
test_XAn.to_csv(args.path_save+'counts_matrix_all_normr_test.tab',sep="\t")

train_XAln.to_csv(args.path_save+'counts_matrix_all_log1p_normr_train.tab',sep="\t")
validation_XAln.to_csv(args.path_save+'counts_matrix_all_log1p_normr_validation.tab',sep="\t")
test_XAln.to_csv(args.path_save+'counts_matrix_all_log1p_normr_test.tab',sep="\t")
