#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Usage: python3 simulate_bulkRS.py --input ../../data/MP/train_valid_test/counts_matrix_pp_test.tab --cell_label_mapping ../../data/MP/train_valid_test/cell_label_mapping.tab --cell_label_oh ../../data/MP/train_valid_test/cell_labels_oh_test.csv --output ../results/simulated_bulkRS/scMP/scMP_test
#Usage: python3 simulate_bulkRS.py --input ../../../RA/other_datasets/HBC_NSR/train_valid_test/counts_matrix_pp_test.tab --cell_label_mapping ../../../RA/other_datasets/HBC_NSR/train_valid_test/cell_label_mapping.tab --cell_label_oh ../../../RA/other_datasets/HBC_NSR/train_valid_test/cell_labels_oh_7celltypes_test.csv --output ../results/simulated_bulkRS/scHB/scHB_test
#Usage: python3 simulate_bulkRS.py --input ../../data/MP/train_valid_test/counts_matrix_pp_test.tab --cell_label_mapping ../../data/MP/train_valid_test/cell_label_mapping.tab --cell_label_oh ../../data/MP/train_valid_test/cell_labels_oh_test.csv --output ../results/simulated_bulkRS/scMP/scMP_test_origcellprop
#Usage: python simulate_bulkRS.py --input ../data/simulated_bulkRS/Ind7.tab --output ../data/simulated_bulkRS/df_Ind7

import pandas as pd
import numpy as np
from numpy import savetxt
import argparse
from itertools import chain
import csv

def sim_bRS(input_file,count,label_map):
    #First simulated dataset 
    #Generate bRS by calculating sums of all columns.. ie. sum up all counts of genes from different cells

    tot_rows = len(input_file)
    cell_types_unique=input_file.Celltype.unique()
    # 1. Total sum per column:
    #df=input_file.drop(['batch','Celltype','Unnamed: 0','Unnamed: 0.1'],axis=1)
    df=input_file.drop(['Celltype','Unnamed: 0'],axis=1)
    #print(df.sum(axis=0))


    df.sum(axis=0).to_csv(args.output+'_Simset'+str(count)+'.csv')

    # 2. Cell-type proportions:
    cell_types_unique=[]
    proportions=[]
    #print(cell_types_unique)
    for key, value in label_map.items():
        #print(element)
        #print(input_file[input_file['Celltype'].str.contains(element)])
        occurences = input_file[input_file['Celltype'].str.contains(value)]
        #print(len(occurences))
        #frac = (len(occurences)*100)/tot_rows
        frac = (len(occurences))/tot_rows
        cell_types_unique.append(value)
        proportions.append(frac)
        #print(frac)

    df_proportions = pd.DataFrame({'Cell_types':cell_types_unique,'Proportion':proportions})
    df_proportions.set_index('Cell_types',inplace=True)
    df_proportions.to_csv(args.output+'_Simset'+str(count)+'_celltype_proportions.csv')

    # 3. Generate matrix of sums for each gene for each cell_type
    #df=input_file.drop(['batch','Unnamed: 0','Unnamed: 0.1'],axis=1)
    df=input_file.drop(['Unnamed: 0'],axis=1)
    data_tot=[]
    data_avg=[]
    for key, value in label_map.items():
        occurences = df[df['Celltype'].str.contains(value)]
        #print(occurences)
        o = occurences.sum(axis=0)
        #print(occurences.sum(axis=0))
        #print(df[df['Celltype'].str.contains(element)].sum())   same as above line
        data_tot.append(o)
        
        #print(len(df[df['Celltype'].str.contains(element)]))    # Get total count of members
        o = o.drop(['Celltype'])
        #print(o/len(df[df['Celltype'].str.contains(element)]))  # Get avg of counts for genes
        avg = o/len(df[df['Celltype'].str.contains(value)])  # Get avg of counts for genes
        data_avg.append(avg)
        #df_celltype_avg.append(occurences.loc[element])

        #df_celltype_avg.append(occurences.loc[element])

    df_celltype_tot=pd.concat(data_tot,axis=1)
    df_celltype_tot.columns=cell_types_unique
    df_celltype_tot.to_csv(args.output+'_Simset'+str(count)+'_matrix_celltype_tot.csv')
    df_celltype_avg=pd.concat(data_avg,axis=1)
    df_celltype_avg.columns=cell_types_unique
    df_celltype_avg.to_csv(args.output+'_Simset'+str(count)+'_matrix_celltype_avg.csv')

    df=input_file.drop(['Unnamed: 0','Celltype'],axis=1)
    return (df.sum(axis=0),df_proportions)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulating Bulk-RNAseq from scRNAseq data')
    parser.add_argument('--input', type=str, default='../data/bulkRS/', help='input file')
    parser.add_argument('--cell_label_mapping', type=str, default='../data/bulkRS/', help='input file')
    parser.add_argument('--cell_label_oh', type=str, default='../data/bulkRS/', help='input file')
    parser.add_argument('--output', type=str, default='../results/decon/', help='output file directory')

    args = parser.parse_args()

    with open(args.input,'r') as f:
        input_file = pd.read_csv(f,sep="\t")
    print(input_file.shape)

    cell_label_mapping = open(args.cell_label_mapping,'r')

    label_map = {}
    for line in cell_label_mapping.readlines():
        print(line)
        values = line.strip('\n').split('\t')
        label_map[int(values[0])]=values[1]
    print(label_map)

    num_cell_types = len(label_map)
    print(num_cell_types)

    cell_label_oh = open(args.cell_label_oh,'r')

    label_all = [] 
    for line in cell_label_oh:
        values = line.strip('\n').split(',')
        maxindex = np.array(values).argmax()
        label_all.append(label_map[maxindex])

    print(label_all)

    #Append the celltypes to input file
    input_file['Celltype']=label_all
    print(input_file)

    #Drop all unclassified cells
    input_file = input_file[~input_file['Celltype'].isin(['Unclassified','unclassified','Unknown','unknown','Other','other','NA','na'])]

    #Generate bulkRS for the whole sample
    #print('Whole sample')
    #sim_bRS(input_file,5);

    #Generate 100 simulations of cell type proportions based on original cell type proportions
    #Original cell type proportions for MP dataset
    #fractions_cell_types = list(np.random.dirichlet((0.0074,0.0249,0.101,0.474,0.115,0.146,0.074,0.022,0.003,0.019,0.004,0.005,0.004),100).flatten())
    #Original cell type proportions for HBC dataset
    #fractions_cell_types = list(np.random.dirichlet(0.252,0.191,0.228,0.071,0.067,0.125,0.063).flatten())

    #Generate 100 randomly simulated datasets
    simulated_bulkRS_sets = []
    proportions_sets = []
    for ct in range(100):
        print('Simulation...')

        #Randomly get a collection of 90% of rows in the input - as starting point for generating simulated datasets
        sample_file = input_file.sample(frac=0.90)
        print(sample_file.shape)

        #Simulate n' random numbers summing up to 1 - in order to indicate the fraction of 'n' cell types
        fractions_cell_types = list(np.random.dirichlet(np.ones(num_cell_types)).flatten())
        print(fractions_cell_types)

        #cell_types_unique = list(set(label_all))
        #print(cell_types_unique)

        #For each cell type, sample those many rows corresponding to the cell type - and repeat this process till we get half the cells
        data = []
        data_rows = 0
        while(True):
            ctr = 0
            for key, value in label_map.items():
                sample_celltype = sample_file[sample_file['Celltype'] == value]
                rows_to_get = fractions_cell_types[ctr] * len(sample_file)

                #if rows_to_get > len(sample_celltype): #get all the rows of that cell type
                #    sample_celltype1 = sample_celltype
                #else:
                sample_celltype1 = sample_celltype.sample(frac=fractions_cell_types[ctr],replace=True)
                #print(sample_celltype1)
                ctr = ctr + 1
                data.append(sample_celltype1)
                data_rows = data_rows + len(sample_celltype1)

            if(data_rows >= len(sample_file)):
                break

        #print(data)

        if(len(data) > 1):
            simulated_file = pd.concat(data)

            print('simulated_file')
            print(simulated_file.shape)
            #print(simulated_file)

            df,proportions = sim_bRS(simulated_file,ct,label_map)
            simulated_bulkRS_sets.append(df)
            proportions_sets.append(proportions)

    print('simulated bulkRS sets')
    #print(pdList)
    #new_df = pd.concat(pdList)
    #new_df.to_csv(args.output+'_all_simulated_bulkRS_sets.csv')

    simulated_bulkRS = pd.concat(simulated_bulkRS_sets,axis=1)
    #print(simulated_bulkRS.T)
    print(simulated_bulkRS.T.shape)
    simulated_bulkRS.T.to_csv(args.output+'_all_simulated_bulkRS_sets.tab',sep="\t")

    all_proportions = pd.concat(proportions_sets,axis=1)
    #all_proportions.drop(['Proportion'],inplace=True)
    #print(all_proportions.T)
    print(all_proportions.T.shape)
    #all_proportions.T.drop(all_proportions.T.columns[0],axis=1,inplace=True)
    all_proportions.T.to_csv(args.output+'_all_proportions_sets.csv',header=False)
    exit()

    # opening the csv file in 'w+' mode 
    file = open(args.output+'_all_simulated_bulkRS_sets.csv', 'w+', newline='') 
  
      # writing the data into the file 
    with file:     
        write = csv.writer(file) 
        write.writerows(simulated_bulkRS_sets) 
