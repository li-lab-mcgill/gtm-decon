#This code merges 2 phi files based on their genes 
#Usage: python3 get_GTEX_tissue_set.py --path_input ../../data/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct_trunc --path_mapping ../../data/GTEX/Pancreas.lst --path_save ../../data/GTEX/GTEx_Pancreas.csv 
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Merge phi sets')
parser.add_argument('--path_phi1', type=str, default='../data/bulkRS/', help='input file',required=True)
parser.add_argument('--path_phi2', type=str, default='../data/bulkRS/', help='input file',required=True)
parser.add_argument('--path_ref_genes', type=str, default='../data/bulkRS/', help='input file',required=True)
parser.add_argument('--path_save', type=str, default='../data/bulkRS/', help='input file',required=True)

args = parser.parse_args()

with open(args.path_phi1,'r') as f:
        df1 = pd.read_csv(f,sep=",")

with open(args.path_phi2,'r') as f:
        df2 = pd.read_csv(f,sep=",")

df_cd = pd.merge(df1, df2, how='outer', on = 'genes')

with open(args.path_ref_genes,'r') as f:
        ref_genes = pd.read_csv(f,header=None)
        #print(genes.values.tolist())
        #print(X['Unnamed:0'].values.tolist())
genes = list(ref_genes.iloc[:,0])
        #Retain only those genes corresponding to genes in the training dataset
pp_genes = [item for item in genes]
#print(pp_genes)

pp = df_cd.loc[df_cd['genes'].isin(pp_genes)]
        #pp = X.loc[X['Unnamed:0'].isin(genes.values.tolist())]
        #print(pp)
        #Rearrange the rows of the dataframe to match the order of the genes in the training dataset - Note, those that are not found in Bulk Set will be set to NaNs
        #Keep only the first occurence for genes that are duplicated
pp.set_index('genes',inplace=True)
pp = pp[~pp.index.duplicated(keep='first')]
        #print(pp)
pp = pp.reindex(pp_genes)
        #print(pp)
pp = pp.reset_index()
pp = pp.fillna(0)

pp.to_csv(args.path_save,index=False)
