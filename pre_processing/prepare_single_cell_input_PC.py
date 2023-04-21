#This code strips down each single cell matrix into a count matrix ie. gene names, cell IDs, counts. All other meta information is removed

#Author: Swapna Seshadri

#Usage: python3 prepare_single_cell_input_PI_Segerstolpe.py --path_input ../../data/PI_Segerstolpe/pancreas_refseq_counts_goodquality_7ind.csv --path_save ../../data/PI_Segerstolpe/

import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc
sc.settings.verbosity = 2
import argparse

parser = argparse.ArgumentParser(description='Preparing scRNAseq input')
parser.add_argument('--path_input', type=str, default='../data/bulkRS/', help='input file')
parser.add_argument('--path_save', type=str, default='../results/decon/', help='output file directory')

args = parser.parse_args()

#pancreas_refseq_counts_goodquality.csv
with open(args.path_input,'r') as f:
    X = pd.read_csv(f,sep=",")

orig_X = X
print('Saving counts matrix...')
X.to_csv(args.path_save+'counts_matrix.tab',sep="\t")

#Generate cell labels to one-hot encoding mapping
#Label map
label_map={'Acinar_cell': 0, 'Ductal_cell_type_1': 1, 'Ductal_cell_type_2': 2, 'Endocrine_cell': 3, 'Endothelial_cell': 4, 'Fibroblast_cell': 5, 'Stellate_cell': 6, 'Macrophage_cell': 7, 'T_cell': 8, 'B_cell': 9}
label_all =  np.array([label_map[i] for i in X['Celltype']])

## Save cell labels as one hot encoding
a = np.array(label_all)
b = np.zeros((a.size, a.max()+1))
b[np.arange(a.size),a] = 1
print('Saving cell labels...')
np.savetxt(args.path_save+'cell_labels_oh.csv',b,delimiter=",")

X.drop(X.columns[0], axis = 1, inplace=True)
X.drop(X.columns[-1], axis = 1, inplace=True)
X.to_csv(args.path_save+'counts_matrix_all.tab',sep="\t")

## Normalize data - for all gene
X = orig_X
X.drop(X.columns[0], axis = 1, inplace=True)
X.drop(X.columns[-1], axis = 1, inplace=True)
genes = X.columns
adata = AnnData(X[genes])
X_norm = sc.pp.normalize_total(adata, target_sum=10000, exclude_highly_expressed=True, inplace=False)['X']
X_norm_r = np.rint(X_norm)
df = pd.DataFrame(data=X_norm_r,index=X.index,columns=X.columns)
df = df.astype(int)
print('Saving all genes normalized counts matrix ..')
df.to_csv(args.path_save+'counts_matrix_all_normr.tab',sep="\t")  

## Log normalize data - for all genes
X_norm_log = sc.pp.log1p(X_norm) #Logarithamize the data
X_norm_log_r = np.rint(X_norm_log)
df = pd.DataFrame(data=X_norm_log_r,index=X.index,columns=X.columns)
df = df.astype(int)
print('Saving all log normalized counts matrix ..')
df.to_csv(args.path_save+'counts_matrix_all_log1p_normr.tab',sep="\t")

#Preprocess the data - simple preprocessing
A = orig_X.T
A = A[A.astype(bool).sum(axis=1) < len(orig_X)*0.8]    #Remove highly expressed genes - ie. found in >80% of cells
A = A[A.astype(bool).sum(axis=1) > 5]             #Remove lowly expressed genes - ie. found in <5 cells
X = A.T
print('Saving preprocessed counts matrix .. ')
X.to_csv(args.path_save+'counts_matrix_pp.tab',sep="\t")

#Normalize the data - for simple preprocessing
genes = X.columns
adata = AnnData(X[genes])
X_norm = sc.pp.normalize_total(adata, target_sum=10000, exclude_highly_expressed=True, inplace=False)['X']
X_norm_r = np.rint(X_norm)
df = pd.DataFrame(data=X_norm_r,index=X.index,columns=X.columns)
df = df.astype(int)
print('Saving preprocessed normalized counts matrix ..')
df.to_csv(args.path_save+'counts_matrix_pp_normr.tab',sep="\t")

#Log normalize the data - after 10000* based normaliation - for simple preprocessing
X_norm_log = sc.pp.log1p(X_norm) #Logarithamize the data
X_norm_log_r = np.rint(X_norm_log)
df = pd.DataFrame(data=X_norm_log_r,index=X.index,columns=X.columns)
df = df.astype(int)
print('Saving preprocessed log normalized counts matrix ..')
df.to_csv(args.path_save+'counts_matrix_pp_log1p_normr.tab',sep="\t")

#Preprocess the data to get only the set of highly variable genes
X = orig_X
#print(X)
X.drop(X.columns[0], axis = 1, inplace=True)
X.drop(X.columns[-1], axis = 1, inplace=True)
genes = X.columns
#print(genes)
adata_raw = AnnData(X[genes])
adata = adata_raw.copy()
#print(adata.X)
sc.pp.normalize_total(adata, exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
print("Highly variable genes: %d"%sum(adata.var.highly_variable))
var_genes_all = adata.var.highly_variable

#Using the overall set of highly variable genes
var_select = adata.var.highly_variable
var_genes = var_select.index[var_select]
print('Mode - overall')
print(len(var_genes))
#print(var_genes)
adata2 = adata_raw[:,var_genes]
#print(adata2)

genes2 = adata2.var.index
df2 = pd.DataFrame(adata2.X, columns=genes2)

#print('genes')
#print(genes2)
#print('df2')
#print(df2)
df2 = df2.astype(int)
print('Saving counts matrix.. highly variable genes... ')
df2.to_csv(args.path_save+'counts_matrix_hvg.tab',sep="\t") #Counts for variable genes

adata3 = adata[:,var_genes]
genes3 = adata3.var.index
X_norm_r = np.rint(adata3.X)
df3 = pd.DataFrame(X_norm_r, columns=genes3)
df3 = df3.astype(int)
print('Saving counts matrix.. highly variable genes log1p... ')
df3.to_csv(args.path_save+'counts_matrix_hvg_log1p_normr.tab',sep="\t") #Normalized  for variable genes
