# GTM-decon: Guided Topic Modeling for Deconvolution of cell types from bulk RNA-seq data

![My Image](main_fig.png)

**GTM-decon model overview:** GTM-decon uses a guided topic modeling approach to model cell-type information from single cell RNAseq data and subsequently use it to deconvolve cell-type proportions from bulk RNAseqRNA-seq data. GTM-decon works by factorizing the gene expression data from the large single cell gene expression datasets (cells-by-genes matrices) into two matrices: a genes-by-topics (i.e., cell types) matrix φ (φ) capturing the importance of gene expression of each gene for each cell type, and a cells-by-topics matrix (θ) capturing the importance of different topics (i.e., cell types) in each cell. The modelling is guided by the prior knowledge of cell type labels in the single cell datasets. Each cell-type is modelled as a topic (or a set of topics), by enforcing topic to cell type assignment by means of the Bayesian topic prior values assigned to the α hyperparameter, which are set to high values in comparison with the rest of the topics for topic(s) corresponding to that cell type (see Methods). This forces the topic inference (θ) to anchor each topic at a cell type and subsequentaly influences the genes-by-topics (φ) inference. The latter as a global parameter is subsequently used to infer the cell-type proportion in a bulk RNA-seq sample, which manifests as averaged gene expression from all cell types in a tissue of interest. Deconvolution of bulk RNA-seq datasets enables their segregation into different constituent cell types, along with the elucidation of cell-type proportions. This in turn can be used as a surrogate to gauge the biological aspects of the sample from which the bulk RNA-seq data was obtained, such as its health (healthy vs. diseased) as well as developmental-state or cell-state.

GTM-decon is a Unix-style command-line tool. You can compile it on a unix machine.

## INSTALLATION:

To install GTM-decon, you will need to first install armadillo (http://arma.sourceforge.net)

Assuming you are in the mixehr-surelda directory, to compile, simply run:
```
make
```

## Input data:
See toy_dataset.tar.gz for format of training data (vased on pancreatic scRNA-seq data)

To run GTM-decon, the following files are required:

1. metaData.txt file - Stores information about the phenotypes ie. genes in our case in the format <typeID, pheID, stateCnt> where typeID indicates a distinct phenotype (in this cases genes - designated 1, pheID corresponds to phenotypeID for each gene, ie. gene1, gene2, ... geneN), and stateCnt indicates number of states for phenotype - corresponds to only one state in this case 
2. trainData.txt file - Stores information about the counts data in the format <cellID, typeID, pheID, stateID, freq> where typeID and pheID correspond to metaData file above, cellID corresponds to the cell IDs, as in 1, 2, 3... M, stateID is 0 based in this case, freq corresponds to value of counts. This information is provided only for those cells and genes where the count is non-zero
3. priorData.txt file - Stores information about prior probabilites for each cell type for each cell ID in the format <cellID, topicId, priorprob> where topicID corresponds to the metaphenotype ID (corresponding to 'N' topics for 'N' cell types in most of the general cases), starting from an index of 0, and the prior probability for the cell type for that cellID (must be >0 and <=1). GTM works by using the one-hot cell type encodings for the cellIDs to assign prior probabilites. The topic corresponding to the cell type is assigned a prior probability of 0.9, and the rest of the cell types are assigned a value corresponding to 0.1/(number of cell types - 1).

Steps involved in the approach:

1. Generation of cell-type specific gene signatures from scRNAseq datasets

   a. Preparing the single cell RNAseq dataset in suitable input formats:

       The key inputs required from single cell RNAseq dataset are:

	- Count matrix (as cell types x genes)
	- Meta data corresponding to "Celltype" for each cellID, and "Batch" information corresponding to the individual / sample from which the cellID was obtained

       For evaluation purposes, we create an artificial bulk RNAseq dataset from the scRNAseq data (by summing up the number of cells belonging to a cell type from the dataset), and derive a matrix of cell-type proportions based on the artificially constructed set.

       A sample python script for generating these required files is as follows: (This should be tweaked according to the dataset)

	python3 bulkRS_from_scRS_PI_Segerstolpe_cellprop.py --input <path to scRNAseq data file in genes x cell types format> --meta <path to meta data file containing cell type information for cellIDs> --output_bRS  <path to which output file containing artificially constructed data should be saved> --output_sc <path to which output file containing filtered scRNAseq data in the correct format should be saved> --sample_order <Sample order in which artificial bulkRS data should be saved> (optional, but preferable) --celltype_order <Order of cell types in the cell proportions file> (optional, but preferable) --output_cp <path to file containing cell type proportions>
