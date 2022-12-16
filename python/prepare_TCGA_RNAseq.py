#This code strips down each single cell matrix into a count matrix ie. gene names, cell IDs, counts. All other meta information is removed

import argparse

parser = argparse.ArgumentParser(description='Preparing TCGA RNAseq data')
parser.add_argument('--path_input', type=str, default='../data/bulkRS/', help='input file',required=True)
parser.add_argument('--path_save', type=str, default='../results/decon/', help='output file directory',required=True)

args = parser.parse_args()

I = open(args.path_input,'r')
O = open(args.path_save,'w')

#Retain only 1) rows with gene names 2) and retain only gene names
#Get the sample ID for each cell ID
for line in I.readlines():
    values = line.strip('\n').split('\t')
    if values[0] != "gene_id" and "?" not in values[0]: #do not consider these rows
            gene_name = values[0].split('|')[0]
            if values[0] == "Hybridization REF":
                print("Unnamed:0",end=",",file=O)
            else:
                print(gene_name,end=",",file=O)
            print(','.join(map(str,values[1:])),file=O)
