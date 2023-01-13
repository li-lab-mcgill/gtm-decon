#This code strips down each single cell matrix into a count matrix ie. gene names, cell IDs, counts. All other meta information is removed
#python3 prepare_bulkRS_hvg.py --path_input ../../data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.TNBC_ER.tab --path_save ../../data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.TNBC_ER
#python3 prepare_bulkRS_hvg.py --path_input ../../data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.uniq.type.2sets.tab --path_save ../../data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.pp.uniq.type.2sets

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
    X = pd.read_csv(f,sep="\t")

orig_X = X

#Preprocess the data to get only the set of highly variable genes
#print(X)
X.drop(X.columns[0], axis = 1, inplace=True)
X.drop(X.columns[-1], axis = 1, inplace=True)
#X.drop(X.columns[-1], axis = 1, inplace=True)
#X.drop(X.columns[-1], axis = 1, inplace=True)
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
print('Saving bulkRS matrix.. highly variable genes... ')
df2.to_csv(args.path_save+'_hvg.tab',sep="\t") #Counts for variable genes

