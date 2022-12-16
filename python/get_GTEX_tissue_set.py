#This code extracts all columns corresponding to a tissue from GTEX
#Usage: python3 get_GTEX_tissue_set.py --path_input ../../data/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct_trunc --path_mapping ../../data/GTEX/Pancreas.lst --path_save ../../data/GTEX/GTEx_Pancreas.csv 
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Extracting GTEX data for a tissue')
parser.add_argument('--path_input', type=str, default='../data/bulkRS/', help='input file',required=True)
parser.add_argument('--path_mapping', type=str, default='../data/bulkRS/', help='input file',required=True)
parser.add_argument('--path_save', type=str, default='../results/decon/', help='output file directory',required=True)

args = parser.parse_args()

with open(args.path_input,'r') as f:
        df = pd.read_csv(f,sep="\t")

mapping = open(args.path_mapping,'r')
map_list = mapping.readlines()
map_list = list(map(str.strip, map_list))
map_list.insert(0,"Description")

tissue_df = df[df.columns.intersection(map_list)]
#print(tissue_df)

tissue_df.to_csv(args.path_save,index=False)
