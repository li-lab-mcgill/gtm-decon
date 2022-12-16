#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This program generates artificial bulkRS from scRS by summing up all the counts for each sample
# Usage: python3 bulkRS_from_scRS.py --input /home/mcb/users/slaksh1/projects/scDECON_CANCER/data/PI_Segerstolpe_scRS/pancreas_refseq_counts_3514sc.txt --meta /home/mcb/users/slaksh1/projects/scDECON_CANCER/data/PI_Segerstolpe_scRS/E-MTAB-5061.sdrf.txt --output_bRS  /home/mcb/users/slaksh1/projects/scDECON_CANCER/data/PI_Segerstolpe_scRS/pancreas_artificial_bulkRS.csv --output_sc  /home/mcb/users/slaksh1/projects/scDECON_CANCER/data/PI_Segerstolpe_scRS/pancreas_refseq_counts_goodquality.csv

import pandas as pd
import argparse
import csv

#Taken from https://thispointer.com/python-how-to-find-keys-by-value-in-dictionary/#:~:text=As%2C%20dict.,in%20a%20separate%20list%20i.e.
def getKeysByValue(dictOfElements, valueToFind):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item  in listOfItems:
        if item[1] == valueToFind:
            listOfKeys.append(item[0])
    return  listOfKeys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generating Bulk-RNAseq from scRNAseq data')
    parser.add_argument('--input', type=str, default='../data/bulkRS/', help='input file')
    parser.add_argument('--meta', type=str, default='../data/bulkRS/', help='input file')
    parser.add_argument('--output_bRS', type=str, default='../results/decon/', help='output file directory')
    parser.add_argument('--output_sc', type=str, default='../results/decon/', help='output file directory')

    args = parser.parse_args()

    with open(args.input,'r') as f:
        scRS_file = pd.read_csv(f,sep="\t")
    #print(scRS_file.shape)

    meta = open(args.meta,'r')

    #Get the sample ID for each cell ID
    sample_mapping_dict = {}
    celltype_mapping_dict = {}
    samples = []
    for line in meta.readlines():
        values = line.strip('\n').split('\t')
        if values[7] != "not applicable":
            sample_mapping_dict[values[0]] = values[0].split('_')[0]
            celltype_mapping_dict[values[0]] = values[7]
            if values[0].split('_')[0] != 'Source Name': 
                samples.append(values[0].split('_')[0])
    uniq_samples = list(set(samples))

    #Get subset of dataframe columns corresponding to each sample - sum up values for each genes for all the sample IDs 

    df1 = pd.DataFrame({"Genes":scRS_file["#samples"]})
    df3 = pd.DataFrame({"Genes":scRS_file["#samples"]})
    #print('genes')
    #print(df1)
    all_sampleIDs=[]
    for sample in uniq_samples:
        sampleIDs = getKeysByValue(sample_mapping_dict,sample)
        all_sampleIDs.extend(sampleIDs)
        #print('sampleIDs')
        #print(sampleIDs)
        #dataframe from scRNAseq corresponding to specific celltype for specific sample - Sum up all gene counts to get bulkRS value
        subset_df = scRS_file[sampleIDs]
        #print('subset_df')
        #print(subset_df)
        df2 = pd.DataFrame({sample:subset_df.sum(axis=1)})
        df1 = df1.join(df2[sample])
        df3 = df3.join(subset_df)
        #print('df2')
        #print(df2)
    #print('df1')
    #print(df1)
    df1.to_csv(args.output_bRS,index=False)

    df4 = df3[0:26177] #to remove spike-in information
    #print(df4)
    celltypes = [celltype_mapping_dict[i] for i in all_sampleIDs]
    all_celltypes = ['Celltype']
    all_celltypes.extend(celltypes)
    all_celltypes_series = pd.Series(all_celltypes, df4.columns)
    df4 = df4.append (all_celltypes_series, ignore_index = True)
    #print(df4)

    df4.T.to_csv(args.output_sc,header=None)
