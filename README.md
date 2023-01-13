# GTM-decon: Guided Topic Modeling for Deconvolution of cell types from bulk RNA-seq data

![My Image](main_fig.png)

**GTM-decon overview:** 
We developed GTM-decon to model cell-type information from scRNA-seq data and subsequently use it to deconvolve bulk RNA-seq data. GTM-decon works by factorizing gene expression data from large single-cell gene expression datasets (cells-by-genes matrices) into two matrices: a genes-by-topics (i.e., cell type) matrix (&Phi;) capturing the importance of gene expression of each gene for each cell type, and a cells-by-topics matrix (&theta;) capturing the importance of different topics (i.e., cell types) in each cell. The modelling is guided by the prior knowledge of cell type labels in the scRNA-seq datasets. Each cell-type is modelled as a topic (or a set of topics). This guides the topic inference (&theta;) to anchor each topic at a specific cell type and subsequently influences the genes-by-topics (&Phi;) inference. The latter as a global parameter is subsequently used to infer the cell-type proportion in a bulk RNA-seq sample, which manifests as averaged gene expression from all cell types in a tissue of interest. Deconvolution of bulk RNA-seq datasets enables their segregation into different constituent cell types, thereby elucidating their cell-type proportions. This in turn can be used as a surrogate to gauge the cognate biological states of the subjects, such as their health status or disease progression.

GTM-decon is a Unix-style command-line tool. You can compile it on a unix machine.

## Installation:

To install GTM-decon, you will need to first install armadillo (http://arma.sourceforge.net)

Assuming you are in the gtm-decon directory, to compile, simply run:
```
make
```

## Preparing data:
See pancreas_dataset.tar.gz for format of training data (based on pancreatic scRNA-seq data). See decon.sh and DEgenes.sh shell scripts for sample scripts.

To run GTM-decon, the following files are required:

a. metaData.txt file - Stores information about the phenotypes ie. genes in the format <typeID, geneID, stateCnt> where typeID indicates a distinct phenotype (in this cases gene expression - designated 1, geneID corresponds to geneID for each gene, ie. gene1, gene2, ... geneN), and stateCnt indicates number of states for phenotype - corresponds to only one state in this case 
b. trainData.txt file - Stores information about the counts data in the format <cellID, typeID, gebeID, stateID, freq> where typeID and geneID correspond to metaData file above, cellID corresponds to the cell IDs, as in 1, 2, 3... M, stateID is 0 based in this case, freq corresponds to value of counts. This information is provided only for those cells and genes where the count is non-zero
c. priorData.txt file - Stores information about prior probabilites for each cell type for each cellID in the format <cellID, topicId, priorprob> where topicID corresponds to the metageneID (corresponding to 'N' topics for 'N' cell types in most of the general cases), starting from an index of 0, and the prior probability for the cell type for that cellID (must be >0 and <=1). GTM works by using the one-hot cell type encodings for the cellIDs to assign prior probabilites. The topic corresponding to the cell type is assigned a prior probability of 0.9, and the rest of the cell types are randomly assigned a value ranging between 0.01-0.1, for a 1-topic-per-celltype. For 2, 3, 4, and 5 topics- per cell-type models, the topics corresponding to the cell type is assigned a prior probability of 0.45, 0.33, 0.225, and 0.18 respectively, whereas the rest of the topics are randomly assigned a value ranging between 0.001-0.01.

These can be generated from single cell RNA-seq (scRNA-seq) count matrices using the following scripts:

Required information from scRNA-seq count matrices:
    Count matrix (as cell types x genes),
    Meta data corresponding to "Celltype" for each cellID, and "Batch" information corresponding to the individual / sample from which the cellID was obtained
```
python3 prepare_single_cell_input_PI_Segerstolpe.py --path_input <path to input file containing scRNAseq data in cells x genes format along with Celltype and Batch columns per cell> --path_save <path to output directory where all the various preprocessing and one-hot encoded cell label files are to be saved>
```
Next, we split scRNAseq data into training, validation and test datasets for training and evaluation purposes in the ratio 70:10:20
```
python3 train_test_validation_split.py --path_input <path to directory containing all preprocessed input files> --path_save <path to output directory where all the split files are to be saved>
```
The validation data / test data / simulated bulk data / real bulk data are all stored in the same format as training data. No prior information is available / provided for these.

From these files, the input files for GTM-decon are generated using these C++ scripts:
For training data:
```
./singleCellInput <path_to_input_files/counts_matrix_pp_train.tab> <path_to_input_files/cell_labels_oh_train.csv> <path_to_output_files/> <number of topics per celltype>
``` 
The output file genes.txt file contains the gene names used in training, which is useful later on to generate bulk RNAseq data based on the same set of genes.

For validation / test data:
```
./singleCellInput_TestData <path_to_input_files/counts_matrix_pp_test.tab> <path_to_output_files/testData.txt>
```

## Training GTM-decon using scRNA-seq data:
```
./gtm-decon -f $scdata -m $scmeta -trp $scprior -k $K -i $niter --inferenceMethod JCVB0 --maxcores 8 --outputIntermediates

Flags are:
	-f: single cell training data file
	-m: meta file
	-trp: prior file
	-i: number of iterations
	-k: number of topics (for 1-topic model, k = number of cell types, for 2-topic model, k = 2 x number of cell types, for 3-topic model, k = 3 x number of cell types etc.)
	-n: inference method (JCVB0)
	--maxcores: maximum number of CPU cores to use
	--outputIntermediates: (whether output intermediate learned parameters for inspection)
```	

## Deconvolution of cell-type proportions in bulk RNA-seq data using trained single cell models:
The key is to ensure that bulk RNAseq data consists of the same genes used in the training data, in the same gene order. The following scripts can be used to transform bulk RNA-seq data accordingly:
```
python3 prepare_bulkRNAseq_input.py --path_input <bulk RNAseq counts input file> --path_save <path to save output files in> --preprocessed_genes <path containing genes used in training set>
```
To convert it into format required for gtm-decon:
```
./bulkRNAseqInput <path to bulkRNAseq input file> <path to bulkRNAseq output file/$bulkdata>
```
To deconvolve bulk RNA-seq data based on trained models,
```
./gtm-decon -m $scmeta -n JCVB0 --newRSSamplesData $bulkdata -k $K \
            --trainedModelPrefix <path to trained gtm-decon files>/trainData_JCVB0_nmar_K$k_iter$niter \
	    --inferNewSampleRSMetagene --inferRSSampleParams_maxiter 100
```
The output *metagene.csv file is a <samples X topics> matrix, containing the deconvolved topic mixture proportions for each of the samples. For a one-topic-per-celltype model, the number of topics corresponds to number of celltypes. For a "n-topic-per-celltype" model, the number of topics corresponds to 'n' X celltypes. The following script is used to generate the <samples X celltypes> matrix from this file, representing the deconvolved cell-type proportions per sample.
```
./Metagene_topics_avg <metagene file from previous step> <number of topics> 
```

## Evaluating deconvolution accuracy:
	
For cases where the cell-type proportions are known (ie. for simulated / artifical datasets), the similarity between the actual cell-type proportions and inferred topic mixture proportions are calculated using 3 metrics: pearson correlation coefficient (PCC) , root mean square deviation (RMSD), and mean absolute deviation (MAD).
```		
./Deconvolution_metrics <path to known cell type proportions file>  <path to inferred meta phenotype file - ie. *_metagene.csv / *_metagene_norm.csv> 
```
To generate simulated dataset from the held-out 20% test dataset of scRNA-seq:
```
python3 simulate_bulkRS.py --input $input_dir//counts_matrix_pp_test.tab --cell_label_mapping $input_dir//cell_label_mapping.tab --cell_label_oh $input_dir//cell_labels_oh_test.csv --output $output_dir/simulated_bulkRS/sim_test_pp
```

## Using GTM-decon for phenotype-guided training of bulk RNA-seq data:
Most of the steps are essentially similar to that used for cell-type-guided training of scRNA-seq data. The same scripts described above can be used for this purpose. The main difference is in the input bulk RNA-seq data matrix, which is sparsified to make it amenable for working with topic models. Refer to Fig8.R script in the R/ directory for sparsification.
	
## Using GTM-decon as a nested guided topic model to identify cell-type specific differentially expressed genes between phenotypes from bulk RNA-seq data:
Here, the phenotype-labels are used as primary-level, and cell-types as secondary-level. To generate the training data to reflect these guides,
```
./singleCellInput_DE $input_dir/counts_matrix_train.tab $input_dir/cell_labels_oh_train.csv $output_dir/ <number of topics per cell-type> <number of cell-types>
```
The approach also requires pre-trained GTM-decon &Phi; matrices, concatenated X <number of phenotypes> to serve as the input &Phi; matrix for training. The key is to transform this matrix to contain the same genes present in the bulk RNAseq data for phenotypes, in the same order. Refer to the script *.sh for an example.

The input cell-type specific pre-trained &Phi; matrices are fine-tuned to reflect changes in the phenotypes using GTM-decon:
```
./gtm-decon --outputIntermediates -f $scdata -m $scmeta -trp $scprior -k $K -i $niter --inferenceMethod JCVB0 --maxcores 10 --presetTopicsPrefix $preset_path

New flag:
	--presetTopicsPrefix $preset_path (path to the concatenated input &Phi; matrix)
```
