# GTM_decon Python Wrapper Documentation

The **GTM_decon** class is a python wrapper built around the GTM-decon C code engine, and is intended to simplify its usage. We can use a GTM_decon object to ingest data for the underlying engine, train our model, and perform cell-type inference.

For examples see the vignettes folder for examples: ```/gtm-decon/vignettes/``` 
- ```tutorial.ipynb``` 
- ```example_loocv_script.py```


## Constructor Arguments
All constructor arguments are optional, but it is highly recommended to set the: **experiment_name**, **n_topics**, and **engine_path** arguments, the rest can generally be ignored.
- **experiment_name** (str) : string to describe the experiment to be run, completely optional, but useful for posterity.
  - defaults to ```empty string```
- **n_topics** (int) : number of topics per cell-type for the model, an integer between 1-5
  - defaults to ```1```
- **engine_path** (str) : *strongly recommended*, path to the compiled C executable ```/gtm-decon/gtm-decon-code/gtm-decon```
  - defaults to ```empty string```
- **genes** (list[str]) : *not recommended*, ingesting single-cell reference will set the genes attribute 
  - we can pass a list of genes to GTM_decon, useful when we already have trained a model before and want to infer proportions on bulk, and didn't use the save_parameters function
  - defaults to ```[]```
- **celltypes** (list[str]) : *not recommended*, ingesting single-cell reference will set the genes attribute
  - ingestion of single-cell reference will set this attribute to a sorted list of unique cell-types in the reference
  - defaults to ```[]```
- **bulk_samples** (list[str]) : *not recommended*, ingesting bulk data will set this attribute
  - list containing the order of the bulk sample names
  - defaults to ```[]```
- **verbose** (bool) : when set to True will print updates during the process
  - defaults to ```True```
- **output_intermediates** (bool) : when set to True will save model parameters at some intermediate iterations instead of just the final model
  - defaults to ```False```
- **override_geneset** (bool) : *not recommended*, when set to True, will tell GTM_decon to use the genes provided in the **genes** argument instead of all of the genes in the single-cell reference.
  - this process makes ingesting single-cell reference data to be much slower, it is highly recommended that we perform feature selection before handing it to GTM_decon, so that the reference only contains the genes we intend to use
  - defaults to ```False```


```python
GTM = GTM_decon(
    experiment_name="README-example",
    n_topics=1,
    engine_path="/path/to/gtm-decon", # C executable
    verbose=True,
    output_intermediates=False,      
)
```


## GTM_decon.pipeline
Main method, will take single-cell reference data and target bulk counts data, ingest them into GTM-decon formatted text files, perform training, and perform cell-type proportion prediction.

Can take input data as either ```pd.DataFrame``` or ```ad.AnnData``` format.
- **bulk_data** (pd.DataFrame | ad.AnnData)
  - ```pd.DataFrame``` : bulk counts DataFrame (*genes x samples*)
    - rows are genes (genes x samples unlike single-cell reference format)
    - columns are bulk sample names
  - ```ad.AnnData``` : bulk counts AnnData
    - counts are either saved as ```AnnData.X``` or in the ```AnnData.layers['Counts']```
    - genes are stored either as ```AnnData.var.index``` or as a column in ```AnnData.var```
- **reference_data** (pd.DataFrame | ad.AnnData)
  - ```pd.DataFrame``` : single-cell counts reference DataFrame (*cells x genes*)
    - rows are cells
    - columns are genes + 1 column named "Celltype" containing cell-type labels
  - ```ad.AnnData``` : single-cell counts reference AnnData
    - counts are either saved as ```AnnData.X``` or in the ```AnnData.layers['Counts']```
    - genes are stored either as ```AnnData.var.index``` or as a column in ```AnnData.var```
- **ref_adata_var_gene_column** (str) : see ```GTM_decon.generate_reference_input```
  - defaults to ```None```
- **ref_adata_obs_celltype** (str) : see ```GTM_decon.generate_reference_input```
  - defaults to ```"Celltype```
- **ref_adata_layers_name** (str) : see ```GTM_decon.generate_reference_input```
  - defaults to ```None```
- **bulk_adata_var_gene_column** (str) : see ```GTM_decon.generate_bulk_input```
  - defaults to ```None```
- **bulk_adata_obs_sample_name** (str) : see ```GTM_decon.generate_bulk_input```
  - defaults to ```None```
- **bulk_adata_layers_name** (str) : see ```GTM_decon.generate_bulk_input```
  - defaults to ```None```
- **n_topics** (int) : number of topics per cell-type to use for pipeline
  - sets ```self.n_topics```, only provide if want to overwrite the current ```self.n_topics```
- **engine_path** (str) : only provide if you wish to change
  - defaults to ```None```
- **n_iter** (int) : number of iterations to train model for
  - defaults to ```5```
- **maxcores** (int) : max number of cpu cores to use
  - defaults to ```10```
- **output_intermediates** (bool) : only provide if you want to overwrite the current ```self.output_intermediates```
  - defaults to ```False```


```python
GTM.pipeline(
    bulk_data = bulk_df, # can provide AnnData here instead
    reference_data = reference_df, # can provide AnnData here instead
    directory = "/path/to/save/results/"
)
```


## GTM_decon.generate_reference_input
Ingests reference single cell data into GTM-decon formatted text files for training. Can provide the data as either a DataFrame or AnnData object (only one of DataFrame / AnnData should be provided).
- **directory** (str) : target directory path to save training txt files
- **reference_DataFrame** (pd.DataFrame) : single-cell counts reference DataFrame (*cells x genes*)
  - rows are cells
  - columns are genes + 1 column named "Celltype" containing cell-type labels
  - defaults to ```None```
- **reference_AnnData** (ad.AnnData) : single-cell counts reference AnnData
  - counts are either saved as ```AnnData.X``` or in the ```AnnData.layers['Counts']```
  - genes are stored either as ```AnnData.var.index``` or as a column in ```AnnData.var```
  - defaults to ```None```
- **adata_var_gene_column** (str) : provided if the gene names we wish to use for AnnData input aren't just the index
  - if using AnnData and is set to None, will use ```AnnData.var.index``` as the gene names for ```GTM_decon.genes```
  - defaults to ```None```
- **adata_obs_celltype** (str) : name of the cell-type information column for AnnData input
  - if using AnnData, the name of the column in ```AnnData.obs``` containing the labeled cell-types
  - defaults to ```"Celltype"```
- **adata_layers_name** (str) : optional layers name 
  - if using AnnData, and the counts aren't stored in ```AnnData.X```, retrieves the information from ```AnnData.layers[adata_layers_name]```
  - defaults to ```None```

```python
# DataFrame example
GTM.generate_reference_input(
    directory = '/path/to/save/txt-files/',
    reference_DataFrame = single_cell_reference_df, # Pandas DataFrame
)

# AnnData example
GTM.generate_reference_input(
    directory = '/path/to/save/txt-files/',
    reference_AnnData = single_cell_reference_adata, # AnnData Object
    adata_var_gene_column = None,
    adata_obs_celltype = 'Celltype',
    adata_layers_name = 'Counts'
)
```

## GTM_decon.generate_bulk_input
Ingests bulk counts data into GTM-decon formatted text file (```bulkData.txt```). Can provide the data as either a DataFrame or an AnnData object.
- **directory** (str) : target directory path to save ```bulkData.txt```
- **bulk_DataFrame** (pd.DataFrame) : bulk counts DataFrame (*genes x samples*)
  - rows are genes (genes x samples unlike single-cell reference format)
  - columns are bulk sample names
  - defaults to ```None```
- **bulk_AnnData** (ad.AnnData) : bulk counts AnnData
  - counts are either saved as ```AnnData.X``` or in the ```AnnData.layers['Counts']```
  - genes are stored either as ```AnnData.var.index``` or as a column in ```AnnData.var```
  - defaults to ```None```
- **adata_var_gene_column** (str) : provided if the gene names we wish to use for AnnData input aren't just the index
  - if using AnnData and is set to None, will use ```AnnData.var.index``` as the gene names
  - defaults to ```None```
- **adata_sample_name** (str) : name of the sample batch names column for AnnData input
  - if using AnnData, the name of the column in ```AnnData.obs``` containing the labeled cell-types, if None will use ```AnnData.obs.index```
  - defaults to ```"None"```
- **adata_layers_name** (str) : optional layers name 
  - if using AnnData, and the counts aren't stored in ```AnnData.X```, retrieves the information from ```AnnData.layers[adata_layers_name]```
  - defaults to ```None```

```python
# DataFrame example
GTM.generate_bulk_input(
    directory = '/path/to/save/txt-file/',
    reference_DataFrame = bulk_df, # Pandas DataFrame
)

# AnnData example
GTM.generate_reference_input(
    directory = '/path/to/save/txt-files/',
    reference_AnnData = single_cell_reference_adata, # AnnData Object
    adata_var_gene_column = None,
    adata_sample_name = 'BatchNames',
    adata_layers_name = 'Counts'
)
```


## GTM_decon.run_training
Trains GTM-decon model on the processed and ingested single-cell reference data.
- **directory** (str) : path to the directory containing GTM-format text files (generated by ```GTM_decon.generate_reference_input```)
- **n_iter** (int) : number of iterations to train for
  - defaults to ```5```
- **maxcores** (int) : number of max cpu cores to use
  - defaults to ```10```
- **output_intermediates** (bool) : flag to save GTM-decon models at intermediate iterations
  - defaults to ```False```
- **inferenceMethod** (str) : don't change this
```python
GTM.run_training(
    directory = '/path/to/directory/with/ingested/reference-data/',
    n_iter = 5,
)
```

## GTM.run_deconvolution
Performs cell-type proportion inference of processed GTM-format bulk data. Should only be run after training has been completed.
- **directory** (str) : path to directory containing trained model
- **bulk_directory** (str) : only used if the processed and ingested bulk data (bulkData.txt) was moved to another directory (*not_recommended*)
  - defaults to ```empty string```
- **maxcores** (int) : number of max cpu cores to use
  - defaults to ```10```
- **inferenceMethod** (str) : don't change this
```python
GTM.run_deconvolution(
    directory = '/path/to/directory/with/trained/GTM-decon/model/'
)
```

## GTM_decon.generate_metagene_normalized
Needs to be run when ```GTM_decon.n_topics > 1``` to convert inferred bulk ```metagene``` file to ```metagene_normalized``` files.
- **directory** (str) : directory where we've already run ```GTM.run_deconvolution```
- **inferenceMethod** (str) : don't touch this

```python
GTM.generate_metagene_normalized(
    directory = "/path/to/deconvolved/metagene/files/"
)
```

## GTM_decon.gather_results
Combines ```metagene``` and ```metagene_normalized``` files from running ```GTM.run_deconvolution``` and ```GTM.generate_metagene_normalized``` (when ```n_topics > 1```). Into a single csv, where the rows are the bulk samples and the columns are the cell-types, with an additional column if ```GTM.output_intermediates == True```.
- **directory** (str) : directory where the ```GTM.run_deconvolution``` output results 
- **output_path** (str) : optional path to save the combined results csv file
  - if not provided will by default save ```gatheredResults.csv``` to the same directory 
  - defaults to ```empty string```
- **inferenceMethod** (str) : don't touch this
```python
GTM.gather_results(
    directory = "/path/to/deconvolved/metagene/files/"
)
```



## GTM_decon.set_engine_path
Setter method for ```GTM_decon.engine_path```
- **GTM_engine_path** (str) : path to the compiled gtm-decon C executable
```python
GTM.set_engine_path(
    GTM_engine_path = '/path/to/gtm-decon' # compiled C code
)
```

## GTM_decon.use_geneset
Setter method for ```GTM_decon.genes``` - *not recommended*

Forces GTM_decon to only use the provided genes when ingesting single-cell reference and bulk. Requires ```GTM_decon.override_geneset == True```

- **genes** (list[str]) : list of genes we intend to use
```python
GTM.use_geneset(
    genes = ['Gene1', 'Gene2', 'Gene3', ...]
)
```

## GTM_decon.load_parameters
Loads parameters of previously run GTM_decon experiment from a JSON file
- **path** (str) : path to directory where the previous experiment saved its parameters
- **json_name** (str) : only used if the JSON file was intentionally renamed, defaults to the empty string, and generally won't be used
  - defaults to ```empty string```
```python
GTM.load_parameters(
    path = '/path/to/directory/containing/GTMParameters.json'
)
```

## GTM_decon.save_parameters
Saves parameters of GTM_decon to a JSON file

- **path** (str) : path to directory to save JSON results
- **json_name** (str) : only used if we want to set the JSON name to something other than the default ```GTMParameters.json```, generally best to not use this
  - defaults to ```empty string```
```python
GTM.save_parameters(
    path = '/path/to/directory/',
)
```


## GTM_decon.generate_reference_input_C_implementation 
*Not recommended.* We included C++ executables that are much faster for ingesting data to text format. However these require the data to be saved in a different format. This requires the counts information to be saved in a tab separated file, and also requires an additional csv file containing the cell-type information in a One-Hot encoded manner.
- **directory** (str) : directory where GTM input text files will be saved
- **path_to_C_file** (str) path to C++ executable ```singleCellInput```
- **input_sc_tab_path** (str): path to single cell reference tab file (*cells x genes*)
  - tab separated counts file, where the rows are cells and columns are genes
- **input_labels_oh_path** (str) : path to one hot encoded cell type labels csv file
  - one hot encodings csv of cell-type labels 1 in the corresponding column if a cell is is of that cell-type, 0 otherwise
- **n_topics** (int) : number of topics per cell-type
  - defaults to ```5```

```python
GTM.generate_reference_input_C_implementation(
    directory = "/path/to/save/txt-files/",
    path_to_C_file = "/path/to/singleCellInput",
    input_sc_tab_path = "/path/to/single_cell_data.tab",
    input_labels_oh_path = "/path/to/cell_labels_one_hot.csv",
    n_topics = 5
)
```

## GTM_decon.generate_bulk_input_C_implementation
*Not recommended*. We included C++ executables that can ingest bulk data to GTM-text format faster. Requires tab separated bulk counts information.
- **directory** (str) : directory where GTM input bulk text data will be saved (```bulkData.txt```)
- **path_to_C_file** (str) : path to C++ executable ```bulkRNAseqInput``` 
- **input_bulk_tab_path** (str) : path to bulk counts tab separated file (*genes x samples*)

```python
GTM.generate_bulk_input_C_implementation(
    directory = "/path/to/save/txt-file/",
    path_to_C_file = "/path/to/bulkRNAseqInput",
    input_bulk_tab_path = "/path/to/bulk_data.tab"
)
```
