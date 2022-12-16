#This code maps the clinical notes of TCGA in the order of the RNAseq data - for proper mapping 
#Usage: python3 get_GTEX_tissue_set.py --path_input ../../data/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct_trunc --path_mapping ../../data/GTEX/Pancreas.lst --path_save out.csv
#Usage: python3 map_tcga_rnaseq_clinical.py --path_clinical ../../data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Clinical_Pick_Tier1.Level_4.2016012800.0.0/BRCA.clin.merged.picked.txt --path_rnaseq ../../data/TCGA/BRCA/gdac.broadinstitute.org_BRCA.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/BRCA.csv --path_save ../../data/TCGA/BRCA/Clinical_notes.csv
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Map clinical notes with RNAseq')
parser.add_argument('--path_clinical', type=str, default='../data/bulkRS/', help='input file',required=True)
parser.add_argument('--path_rnaseq', type=str, default='../data/bulkRS/', help='input file',required=True)
parser.add_argument('--path_save', type=str, default='../results/decon/', help='output file directory',required=True)

args = parser.parse_args()

with open(args.path_clinical,'r') as f:
    Cdf = pd.read_csv(f,sep="\t")

with open(args.path_rnaseq,'r') as f:
    Rdf = pd.read_csv(f,sep=",")

rdf_cols = list(Rdf.columns)
rdf_cols.pop(0)
#print(rdf_cols)

matching_names = []
for val in rdf_cols:
    values = val.split('-')
    #print(values)
    val1 = "-"
    name = (values[0],values[1],values[2])
    #print(name)
    val2 = val1.join(name)
    matching_names.append(val2.lower())

#print(matching_names)

tCdf = Cdf.set_index('Hybridization REF').T

new_Cdf = tCdf.reindex(matching_names)

#print(new_Cdf)
new_Cdf.to_csv(args.path_save)
