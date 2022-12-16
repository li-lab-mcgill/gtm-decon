#This code converts a bulk RNAseq matrix consisting of genes x RNAseq data as rows vs columns into a RNAseq x genes count matrix, delimited by tabs.. to be used as "test" or "new patient data" input with MixEHR-sureLDA-scETM

#python3 prepare_bulkRNAseq_input.py --path_input ~/projects/RA/main_datasets/MH.SARDs.RNASeq.COUNT.csv --path_save ~/projects/RA/main_datasets/MH.SARDs.RNASeq.DC --preprocessed_genes ../scHB_labels_7celltypes/scHB_pp_genes.txt
#python3 prepare_bulkRNAseq_input.py --path_input ~/projects/RA/main_datasets/MH.SARDs.RNASeq.COUNT.csv --path_save ~/projects/RA/main_datasets/MH.SARDs.RNASeq.DC.hvg --preprocessed_genes ../results/scHB/scHB_tvt_hvg_7celltypes/scHB_pp_genes.txt
#python3 prepare_bulkRNAseq_input.py --path_input ~/projects/RA/main_datasets/MH.SARDs.RNASeq.COUNT.csv --path_save ~/projects/RA/main_datasets/MH.SARDs.RNASeq.DC --preprocessed_genes ../results/scHB/scHB_tvt_pp_7celltypes/scHB_pp_genes.txt

import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc
import argparse
sc.settings.verbosity = 2

parser = argparse.ArgumentParser(description='Prepare bulk RNAseq count matrix')
parser.add_argument('--path_input', type=str, default='../data/BN/', help='path to input files')
parser.add_argument('--preprocessed_genes', type=str, default='../data/BN/', help='genes from training set only to be included')
parser.add_argument('--path_save', type=str, default='../data/BN/', help='path to save results')
parser.add_argument('--fpkm', type=int, default=0, help='Is the data in FPKM or counts?')

args = parser.parse_args()

with open(args.path_input,'r') as f:
    X = pd.read_csv(f,sep=",")

#Converts all float values if any to int - for bulk RNAseq data - for proper processing by MixEHR-sureLDA
float_col = X.select_dtypes(include=['float64'])
for col in float_col.columns.values:
    X[col] = X[col].astype('int64')

#print(X)
with open(args.preprocessed_genes,'r') as f:
    genes = pd.read_csv(f,header=None)
#print(genes.values.tolist())
#print(X['Unnamed:0'].values.tolist())

#Retain only those genes corresponding to genes in the training dataset
pp_genes = [item for sublist in genes.values.tolist() for item in sublist]
#print(pp_genes)

pp = X.loc[X['Unnamed:0'].isin(pp_genes)]
#pp = X.loc[X['Unnamed:0'].isin(genes.values.tolist())]
#print(pp)
#Rearrange the rows of the dataframe to match the order of the genes in the training dataset - Note, those that are not found in Bulk Set will be set to NaNs
#Keep only the first occurence for genes that are duplicated
pp.set_index('Unnamed:0',inplace=True)
pp = pp[~pp.index.duplicated(keep='first')]
#print(pp)
pp = pp.reindex(pp_genes)
#print(pp)
pp = pp.reset_index()
pp = pp.fillna(0)


#Transpose the data - convert it into "RNAseq_sets" vs "genes" format
if args.fpkm == 1:
    #A = pp.T
    #print(pp)
    s = pp.select_dtypes(include=[np.number])*100
    pp[s.columns] = s
    #print(pp)
    C = pp.round(decimals=0)
    #print(C)
    A = C.T
    #print(A)
    #A = A.astype(int)
    print('Saving FPKM files for bulk RNAseq data.. based on training genes')
    A.to_csv(args.path_save + ".tab",sep="\t",header=0) 
else:
    A = pp.T
    print('Saving counts matrix for bulk RNAseq data.. based on training genes')
    #A = A.astype(int)
    A.to_csv(args.path_save + ".tab",sep="\t",header=0)

    #Normalize the data - for simple preprocessing
    X = pp.T
    X.columns = X.iloc[0]
    X = X[1:]
    #print(X)
    genes = X.columns
    adata = AnnData(X[genes])
    X_norm = sc.pp.normalize_total(adata, target_sum=10000, exclude_highly_expressed=True, inplace=False)['X']
    #X_norm = sc.pp.normalize_total(adata, target_sum=10000, inplace=False)['X']
    X_norm_r = np.rint(X_norm)
    df = pd.DataFrame(data=X_norm_r,index=X.index,columns=X.columns)
    #df = df.astype(int)
    print('Saving preprocessed normalized counts matrix for bulk RNAseq data .. based on training genes')
    df.to_csv(args.path_save+'_normr.tab',sep="\t")

    #Log normalize the data - for simple preprocessing
    genes = X.columns
    adata = AnnData(X[genes])
    sc.pp.log1p(adata) #Logarithamize the data
    X_norm = sc.pp.normalize_total(adata, target_sum=10000, exclude_highly_expressed=True, inplace=False)['X']
    #X_norm = sc.pp.normalize_total(adata, target_sum=10000, inplace=False)['X']
    X_norm_r = np.rint(X_norm)
    df = pd.DataFrame(data=X_norm_r,index=X.index,columns=X.columns)
    #df = df.astype(int)
    print('Saving preprocessed log normalized counts matrix for bulk RNAseq data.. based on training genes')
    df.to_csv(args.path_save+'_log1p_normr.tab',sep="\t")
