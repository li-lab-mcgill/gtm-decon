
#This code converts a phi matrix - to have only lines corresponding to a gene list ie. from bulkRS file

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

args = parser.parse_args()

with open(args.path_input,'r') as f:
    X = pd.read_csv(f,sep=",")

#print(X)
with open(args.preprocessed_genes,'r') as f:
    genes = pd.read_csv(f,header=None)
#print(genes.values.tolist())
#print(X['Unnamed:0'].values.tolist())

#Retain only those genes corresponding to genes in the training dataset
pp_genes = [item for sublist in genes.values.tolist() for item in sublist]
#print(pp_genes)

pp = X.loc[X[X.columns[0]].isin(pp_genes)]
#pp = X.loc[X['Unnamed:0'].isin(genes.values.tolist())]
#print(pp)
#Rearrange the rows of the dataframe to match the order of the genes in the training dataset - Note, those that are not found in Bulk Set will be set to NaNs
#Keep only the first occurence for genes that are duplicated
pp.set_index(X.columns[0],inplace=True)
pp = pp[~pp.index.duplicated(keep='first')]
#print(pp)
pp = pp.reindex(pp_genes)
#print(pp)
pp = pp.reset_index()
pp = pp.fillna(0)

##Adding a column of 1s on the left.. to match phi format
pp.insert(0,'datatype','1')

print('Saving phi matrix.. based on training genes')
pp.index = np.arange(1, len(pp)+1)
pp.to_csv(args.path_save + ".csv",sep=",",header=0)
