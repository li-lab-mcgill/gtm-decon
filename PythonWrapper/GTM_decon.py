#!python3


import json
import pandas as pd
import warnings
import random
import subprocess
import os
import re
import anndata

RAND_UNIF_LOWER_BOUND_1TOPIC = 0.01
RAND_UNIF_UPPER_BOUND_1TOPIC = 0.1
RAND_UNIF_LOWER_BOUND_NTOPIC = 0.001
RAND_UNIF_UPPER_BOUND_NTOPIC = 0.01

REFERENCE_FILE = "trainData.txt"
META_FILE = "meta.txt"
GENES_FILE = "genes.txt"
PRIOR_FILE = 'priorData.txt'
BULKDATA_FILE = "bulkData.txt"
GATHERED_RESULTS_FILE = "gatheredResults.csv"
PARAMETERS_FILE = "GTMParameters.json"

class GTM_decon:
    '''
    Wrapper class for GTM-decon engine.

    Attributes:
        experiment_name : str
        n_topics : int
        engine_path : str
            - path to GTM-decon engine
        genes : list[str]
        celltypes : list[str]
            - expected to be sorted alphabetically
        bulk_samples : list[str]
            - names of the bulk samples in original order
        verbose : bool
        output_intermediates : bool
        override_geneset : bool
            - flag to use self.genes instead of the genes in the reference
            - not recommended, better to preprocess genes in DataFrame or AnnData before passing to GTM-decon
    '''

    prior_dict = {
            1: 0.9,
            2: 0.45,
            3: 0.33,
            4: 0.22,
            5: 0.18,
        }

    def __init__(self,
            experiment_name="",
            n_topics=1,
            engine_path=None,
            genes=[],
            celltypes=[],
            bulk_samples=[],
            verbose=True,
            output_intermediates=False,      
            override_geneset=False,   
        ):
        # TODO check load/save param order
        self.experiment_name = experiment_name
        self.n_topics = n_topics
        self.engine_path = engine_path
        self.genes = genes
        self.celltypes = celltypes
        self.bulk_samples = bulk_samples
        self.verbose = verbose
        self.output_intermediates = output_intermediates
        self.override_geneset = override_geneset

        # warning to ensure engine_path is set before running training/decon
        if self.engine_path is None:
            warnings.warn('Engine path not set.')
    

    def set_engine_path(self, GTM_engine_path=None):
        ''' (str -> None)
        Setter method for engine path.

        Args: 
            GTM_engine_path : str
                - file path to GTM-decon engine
        
        Returns:
            None
        '''

        self.engine_path = GTM_engine_path

        assert os.path.exists(self.engine_path), f"Engine path: {self.engine_path} does not exist!"
    

    def __str__(self):
        s = 'GTM-decon wrapper object with attributes:\n'
        for attr_name, attr_value in self.__dict__.items():
            if isinstance(attr_value, list):
                if len(attr_value) > 5:
                    s += f'  - {attr_name} len: {len(attr_value)}: {attr_value[:5]}...\n'
                else:
                    s += f'  - {attr_name}: {attr_value}\n'
            else:
                s += f'  - {attr_name}: {attr_value}\n'

        return s


    def use_geneset(self, genes):
        '''
        Setter method for self.genes. Only used when override_geneset flag is True. Not recommended.

        Args:
            genes : list[str]
        
        Returns:
            None
        '''

        assert self.override_geneset, "Cannot use geneset when self.override_geneset is False"
        
        if self.genes != []:
            warnings.warn("Overriding existing genes list!")

        self.genes = genes


    def load_parameters(self, path, json_name=""):
        '''
        Loads GTM-decon parameters from GTMParameters.json file that was previously saved. Defaults to finding "GTMParameters.json" in path, can optionally provide filename if it is different.

        Args:
            path : str
            json_name : str [optional]
                - only use this if the parameters were previously saved to a file not named "GTMParameters.json" 
                - not recommended
        
        Returns:
            None
        '''

        if json_name: # if user changed parameters json file name
            path = os.path.join(path, json_name)
        else: # default behavior
            path = os.path.join(path, PARAMETERS_FILE)

        try:
            with open(path, 'r') as f:
                attrs = json.load(f)
        except FileNotFoundError:
            raise FileNotFoundError(f'Unable to find to find {path}')

        if self.verbose:
            print('Loading attributes ...')

        for k, v in attrs.items():
            setattr(self, k, v)

            if self.verbose:
                if isinstance(v, list):
                    p = str(v[:min(len(v), 5)])
                    if len(v) > 5:
                        p += " ..."

                    print(f'Setting self.{k} to {p}')
                else:
                    print(f'Setting self.{k} to {v}')

        return None


    def save_parameters(self, path, json_name=""):
        '''
        Saves GTM-decon parameters to json file. Defaults to "GTMParameters.json" in path directory, can optionally provide filename if it is different.

        Args:
            path : str
            json_name : str [optional]
                - only use this if the parameters were previously saved to a file not named "GTMParameters.json" 
                - not recommended
        
        Returns:
            None
        '''

        if json_name:
            path = os.path.join(path, json_name)
        else:
            path = os.path.join(path, PARAMETERS_FILE)

        # keep synced with changes to init
        attrs = dict(
            experiment_name = self.experiment_name,
            n_topics = self.n_topics,
            engine_path = self.engine_path,
            genes = self.genes,
            celltypes = self.celltypes,
            bulk_samples = self.bulk_samples,
            verbose = self.verbose,
            output_intermediates = self.output_intermediates,
            override_geneset = self.override_geneset,
        )

        if self.verbose:
            print(f'Saving attributes {attrs.keys()} to path: {path}')

        try:
            with open(path, 'w') as f:
                json.dump(attrs, f)

        # not fatal if error is hit, but is bad and should warn user
        except FileNotFoundError:
            warnings.warn(f'Unable to write to {path}')
        except:
            warnings.warn('Unable to save GTMWrapper parameters to json.')
            

        if self.verbose:
            print(f'Successfully saved GTMWrapper parameters to {path}')
    

    def run_training(self, directory, n_iter=5, maxcores=10, output_intermediates=False, inferenceMethod="JCVB0"):
        '''
        Calls GTM-decon engine for training.

        Args:
            directory : str
                - path to directory containing text formatted files (trainData.txt, priorData.txt, meta.txt, genes.txt)
            n_iter : int [optional]
                - number of iterations to train for, default 5
            maxcores : int [optional]
                - max number of cores GTM-decon engine will have access to, default 10
            output_intermediates : bool [optional]
                - if True will save model parameters for intermediate steps, default False
                - if EITHER self.output_intermediates or this argument is True, will output intermediate files
            inferenceMethod : str [optional]
                - don't change this

        Returns:
            None
        '''

        # check that engine_path is set and correct
        assert self.engine_path is not None, f"The GTM Engine Path has not been set."
        assert os.path.exists(self.engine_path), f"The GTM Engine Path {self.engine_path} does not exist."
        assert len(self.celltypes) > 0, "self.celltypes has not been set."
        # assertion checks for txt input files
        assert os.path.exists(os.path.join(directory, REFERENCE_FILE)), f'Reference {REFERENCE_FILE} was not found!'
        assert os.path.exists(os.path.join(directory, META_FILE)), f'Reference {META_FILE} was not found!'
        assert os.path.exists(os.path.join(directory, PRIOR_FILE)), f'Reference {PRIOR_FILE} was not found!'

        K = len(self.celltypes) * self.n_topics

        if self.verbose:
            print('Running Inference Engine ...')

        subprocess_commands = [
                self.engine_path,
                "-f", os.path.join(directory, REFERENCE_FILE), # directory/REFERENCE_FILE}
                "-m", os.path.join(directory, META_FILE), # directory/META_FILE}
                "-trp", os.path.join(directory, PRIOR_FILE), # directory/PRIOR_FILE}
                "-k", f"{K}",
                "-i", f"{n_iter}",
                "--inferenceMethod", inferenceMethod,
                "--maxcores", f"{maxcores}",
            ]

        if output_intermediates or self.output_intermediates:
            subprocess_commands.append("--outputIntermediates")
        
        if self.verbose:
            stdout = None
        else:
            stdout = subprocess.DEVNULL

        subprocess.run(subprocess_commands, stdout=stdout)

        if self.verbose:
            print('Completed Running Inference Engine.')


    def run_deconvolution(self, directory, data_directory="", maxcores=10, inferenceMethod="JCVB0"):
        '''
        Runs GTM-decon engine to deconvolve bulk data. Assumes engine has already been trained.

        Args:
            directory : str
                - path to directory where GTM-decon was trained (the same directory as the run_training was pointed to)
            data_directory : str [optional]
                - shouldn't be used unless bulkData.txt was moved somewhere else
                - not recommended to change
            max_cores : int [optional]
                - max number of cores GTM-decon engine will use during deconvolution
            inferenceMethod : str [optional]
                - Do not change this 

        Returns:
            None
        '''

        # check that engine_path is set and correct
        assert self.engine_path is not None, f"The GTM Engine Path has not been set."
        assert os.path.exists(self.engine_path), f"The GTM Engine Path {self.engine_path} does not exist."
        assert len(self.celltypes) > 0, "self.celltypes has not been set."
        # assertion checks for txt input files
        assert os.path.exists(os.path.join(directory, META_FILE)), f'Reference {META_FILE} was not found!'

        if data_directory:
            assert os.path.exists(os.path.join(data_directory, BULKDATA_FILE)), f'Target Bulk {BULKDATA_FILE} was not found!'
        else:
            assert os.path.exists(os.path.join(directory, BULKDATA_FILE)), f'Target Bulk {BULKDATA_FILE} was not found!'


        K = self.n_topics * len(self.celltypes)
        # gather models in directory
        model_list = os.listdir(directory)
        model_list = [re.search(f"{REFERENCE_FILE[:-4]}_{inferenceMethod}_nmar_K{K}_iter\d+_metagene.csv", m) for m in model_list]
        model_list = sorted(set([m.group(0).replace("_metagene.csv", "") for m in model_list if m]))

        # generally data_directory should be the same as directory
        if not data_directory:
            data_directory = directory

        if self.verbose:
            stdout = None
        else:
            stdout = subprocess.DEVNULL

        # deconvolve artificial bulk data
        for model in model_list:
            subprocess.run(
                [
                    self.engine_path,
                    "-m", os.path.join(directory, META_FILE), # {directory/META_FILE}
                    "-n", inferenceMethod,
                    "--newRSSamplesData", os.path.join(data_directory, BULKDATA_FILE), # data_directory/BULKDATA_FILE}
                    "--trainedModelPrefix", os.path.join(directory, model), # directory/model
                    "-k", f"{K}",
                    "--inferNewSampleRSMetagene",
                    "--inferPatParams_maxiter", "100", # ? 
                    "--maxcores", f"{maxcores}"
                ],
                stdout=stdout
            )


    def generate_reference_input(self,
        directory, # required
        reference_DataFrame=None,
        reference_AnnData=None,
        adata_var_gene_column=None,
        adata_obs_celltype='Celltype',
        adata_layers_name=None,
        ):

        '''
        Generates input txt files for GTM-decon. Can take either a dataframe or anndata as input.

        Args:
            directory : str
                - output directory where all GTM formatted txt files will be saved
            reference_DataFrame : pd.DataFrame [optional]
                - cells x genes
                - columns should be genes, with one column for Celltype labels
                - no other columns should be in the DataFrame
            reference_AnnData : ad.AnnData [optional]
                - obs should contain celltype information
                - var must have gene information either as index or in the columns
            adata_var_gene_columns : str [optional]
                - when provided will use AnnData.var[adata_var_gene_columns]
                - defaults to None, when not provided will use var.index as genes
            adata_obs_celltype : str [optional]
                - name of the column containing celltype information in AnnData.obs
            adata_layers_name : str [optional]
                - if there is a specific layer that stores the information instead of AnnData.X
                - ie if counts is stored in AnnData.layers['Counts']
                - defaults to None, uses AnnData.X as if they are counts
        
        Returns:
            None
        '''

        assert (reference_DataFrame is None) ^ (reference_AnnData is None) , "GTMWrapper.generate_reference_input expects only one format of input at a time."

        # TODO handle duplicate genes
        if reference_DataFrame is not None:
            if 'Celltype' not in reference_DataFrame.columns:
                raise Exception("Celltype not found in reference DataFrame columns.")


            celltype_labels = list(reference_DataFrame['Celltype'])
            celltypes = sorted(set(celltype_labels)) # alphabetical
            reference_DataFrame = reference_DataFrame.drop(columns = ['Celltype'])
            genes = list(reference_DataFrame.columns)

        elif reference_AnnData is not None:
            # get genes list from adata.var
            if adata_var_gene_column is not None:
                assert adata_var_gene_column in reference_AnnData.var.columns, 'The provided var column name for genes does not exist!'
                genes = list(reference_AnnData.var[adata_var_gene_column])
            else:
                genes = list(reference_AnnData.var.index)


            # get celltypes list from adata.obs
            assert adata_obs_celltype in reference_AnnData.obs.columns, 'The provided obs celltype column name does not exist!'
            celltype_labels = list(reference_AnnData.obs[adata_obs_celltype])
            celltypes = sorted(set(celltype_labels))

            if adata_layers_name is not None:
                assert adata_layers_name in reference_AnnData.layers, f"The provided layer was not found in {reference_AnnData}"
                reference_AnnData.X = reference_AnnData.layers[adata_layers_name]
            
        else: # shouldn't ever get here
            raise Exception('Need to provide either a reference AnnData or reference DataFrame')

        # set self.genes to be genes in reference if no overriding geneset is provided
        if not self.override_geneset:
            if self.genes != genes and len(self.genes) != 0:
                warnings.warn('Overwriting previous self.genes list!')
            self.genes = genes

        assert len(self.genes) > 0, "self.genes is empty!"
        assert len(set(self.genes).intersection(set(genes))) > 0, "No genes matching between provided gene set and reference!"

        
        if self.celltypes != celltypes and len(self.celltypes) != 0:
            warnings.warn('Overwriting previous self.celltypes list!')
            
        self.celltypes = celltypes

        # genes.txt
        if self.verbose:
            print(f'Saving genes file to {os.path.join(directory, GENES_FILE)} ...')

        with open(os.path.join(directory, GENES_FILE), 'w') as f:
            for g in self.genes:
                print(g, file=f)

        if self.verbose:
            print(f'Successfully wrote genes file to {os.path.join(directory, GENES_FILE)}')


        # meta.txt
        if self.verbose:
            print(f'Saving meta file to {os.path.join(directory, META_FILE)} ...')

        with open(os.path.join(directory, META_FILE), 'w') as f:
            for i in range(1, len(self.genes) + 1):
                print(f"1 {i} 0", file=f)
        
        if self.verbose:
            print(f'Successfully wrote meta file to {os.path.join(directory, META_FILE)}')

        
        # trainData.txt
        if self.verbose:
            print(f'Saving training file to {os.path.join(directory, REFERENCE_FILE)} ...')


        if reference_DataFrame is not None:
            with open(os.path.join(directory, REFERENCE_FILE), 'w') as f:
                for cell_index in range(len(reference_DataFrame)):
                    expr = reference_DataFrame.iloc[cell_index]

                    for gene_index, gene in enumerate(self.genes):
                        # Skip genes not in self.genes
                        if gene not in expr:
                            continue

                        count = expr[gene]
                        if count > 0:
                            print(f"{cell_index + 1} 1 {gene_index + 1} 0 {count}", file=f)
                            # format is: cellindex 1 geneindex 0 count
        
        elif reference_AnnData is not None:
            with open(os.path.join(directory, REFERENCE_FILE), 'w') as f:
                for cell_index in range(reference_AnnData.shape[0]):
                    expr = reference_AnnData.X[cell_index, :]

                    if self.override_geneset: # need to check for each gene in reference_AnnData if it belongs
                        for gene_index, gene in enumerate(self.genes):
                            if gene not in genes:
                                continue
                            else:
                                count = expr[genes.index(gene)]

                            if count > 0:
                                print(f"{cell_index + 1} 1 {gene_index + 1} 0 {count}", file=f)

                    else: # when not using an overriding geneset, can just enumerate self.genes
                        for gene_index in range(len(expr)):
                            count = expr[gene_index]

                            if count > 0:
                                print(f"{cell_index + 1} 1 {gene_index + 1} 0 {count}", file=f)

        else: # shouldn't ever get here
            raise Exception('Need to provide either a reference AnnData or reference DataFrame')
        

        if self.verbose:
            print(f'Successfully wrote training file to {os.path.join(directory, REFERENCE_FILE)}')


        # priorData.txt
        if self.verbose:
            print(f'Saving prior file to {os.path.join(directory, PRIOR_FILE)} ...')

        prior_val = self.prior_dict[self.n_topics]
        with open(os.path.join(directory, PRIOR_FILE), 'w') as f:
            for cell_index, cell in enumerate(celltype_labels):
                for celltype_index, cell_ in enumerate(celltypes):
                    if self.n_topics == 1:
                        if cell == cell_:
                            print(f"{cell_index + 1} {celltype_index} {prior_val}", file=f)
                        else:
                            r = random.uniform(RAND_UNIF_LOWER_BOUND_1TOPIC, RAND_UNIF_UPPER_BOUND_1TOPIC)
                            print(f"{cell_index + 1} {celltype_index} {round(r, 7)}", file=f)
                    else:
                        r = random.uniform(RAND_UNIF_LOWER_BOUND_NTOPIC, RAND_UNIF_UPPER_BOUND_NTOPIC)
                        for j in range(self.n_topics):
                            if cell == cell_:
                                print(f"{cell_index + 1} {celltype_index * self.n_topics + j} {prior_val}", file=f)
                            else:
                                print(f"{cell_index + 1} {celltype_index * self.n_topics + j} {round(r, 7)}", file=f)

        if self.verbose:
            print(f'Successfully wrote prior file to {os.path.join(directory, PRIOR_FILE)}')
    

    def generate_bulk_input(self, 
        directory, 
        bulk_DataFrame=None, 
        bulk_AnnData=None,
        adata_var_gene_column=None,
        adata_obs_sample_name=None,
        adata_layers_name=None,
        ):
        '''
        Args:
            directory : str
                - output directory where GTM formatted txt bulk file will be saved
            bulk_DataFrame : pd.DataFrame [optional]
                - genes x samples
                - index should be genes
                - column name should represent the batches
            bulk_AnnData : ad.AnnData [optional]
                - obs should contain batch information
                - var must have gene information either as index or in the columns
            adata_var_gene_columns : str [optional]
                - when provided will use AnnData.var[adata_var_gene_columns]
                - defaults to None, when not provided will use var.index as genes
            adata_obs_sample_name : str [optional]
                - name of the column containing batch sample names in AnnData.obs
                - defaults to None, when not provided will use AnnData.obs.index as sample names
            adata_layers_name : str [optional]
                - if there is a specific layer that stores the information instead of AnnData.X
                - ie if counts is stored in AnnData.layers['Counts']
                - defaults to None, uses AnnData.X as if they are counts

        Returns:
            None
        '''

        assert (bulk_DataFrame is None) ^ (bulk_AnnData is None) , "GTMWrapper.generate_bulk_input expects only one format of input at a time."
        assert len(self.genes) != 0, "GTMWrapper self.genes has not been intialized!"
            
        
        
        if bulk_DataFrame is not None:
            # force converts float values to ints
            bulk_DataFrame = bulk_DataFrame.astype('int')
            sample_names = list(bulk_DataFrame.columns)

            # handle duplicate indices
            if len(bulk_DataFrame.index) != len(bulk_DataFrame.index.unique()):
                bulk_DataFrame = bulk_DataFrame[~bulk_DataFrame.index.duplicated(keep='first')]
                warnings.warn('Duplicate genes found in bulk data indices! Dropping duplicates, keeping first!')
                    
        
            if len(set(self.genes).intersection(set(bulk_DataFrame.index))) == 0:
                raise Exception('No genes found overlapping between Reference and Bulk data.')
        
        elif bulk_AnnData is not None:
            # if a layers column is provided, use adata.layers counts instead of X
            if adata_layers_name is not None:
                assert adata_layers_name in bulk_AnnData.layers, f'adata_layers_name was provided but not found in {bulk_AnnData}'
                bulk_AnnData.X = bulk_AnnData.layers[adata_layers_name]

            bulk_AnnData.X = bulk_AnnData.X.astype('int')

            if adata_obs_sample_name is not None:
                if adata_obs_sample_name not in bulk_AnnData.obs.columns:
                    warnings.warn('Using bulk_AnnData.obs.index as sample names')
                    sample_names = list(bulk_AnnData.obs.index)
                else:
                    sample_names = list(bulk_AnnData.obs[adata_obs_sample_name])
            else:
                # default uses obs.index as sample names
                sample_names = list(bulk_AnnData.obs.index) 
            
            # WARNING: duplicate genes are not handled in adata format, assumes they've been cleaned already


        # Stores sample names in self.bulk_samples
        self.bulk_samples = sample_names


        if self.verbose:
            print(f"Saving bulk text format file to {directory} ...")

        if bulk_DataFrame is not None:
            with open(os.path.join(directory, BULKDATA_FILE), 'w') as f:
                # loop over samples
                for sample_index, sample in enumerate(bulk_DataFrame.columns):
                    # loop over genes
                    for gene_index, gene in enumerate(self.genes):
                        if gene in bulk_DataFrame.index:
                            count = bulk_DataFrame.loc[gene][sample]
                            if count > 0:
                                print(f'{sample_index + 1} 1 {gene_index + 1} 0 {count}', file=f)


        elif bulk_AnnData is not None:
            if adata_var_gene_column is not None and adata_var_gene_column in bulk_AnnData.var.columns:
                bulk_genes = bulk_AnnData.var[adata_var_gene_column]
            else:
                bulk_genes = bulk_AnnData.var.index

            with open(os.path.join(directory, BULKDATA_FILE), 'w') as f:
                # loop over samples
                for sample_index, sample in enumerate(sample_names): 
                    # loop over genes
                    for gene_index, gene in enumerate(self.genes): 
                        if gene in bulk_genes:
                            count = bulk_AnnData[sample, gene].X[0][0]
                            if count > 0:
                                print(f'{sample_index + 1} 1 {gene_index + 1} 0 {count}', file=f)
        
        if self.verbose:
            print(f"Successfully wrote bulk file to {os.path.join(directory, BULKDATA_FILE)}")

 
    def generate_reference_input_C_implementation(self, directory, path_to_C_file, input_sc_tab_path, input_labels_oh_path, n_topics=5):
        '''
        If data is pre-formated in tab separated format with OH cell labels file, can use much quicker C implementation to generate txt input files.

        Args:
            directory : str
                - directory where GTM input text files will be saved
            path_to_C_file : str
                - path to C executable
            input_sc_tab_path : str
                - path to single cell reference tab file
            input_labels_oh_path : str
                - path to one hot encoded cell type labels csv file
            n_topics : int [optional]
                - # of topics per celltype
                - defaults to 5

        Returns:
            None
        '''

        assert os.path.exists(path_to_C_file), "Path to singleCellInput C file does not exist."
        assert os.path.isdir(directory), "Path to output directory does not exist."
        assert os.path.exists(input_sc_tab_path), "Path to input single cell tab data does not exist."
        assert os.path.exists(input_labels_oh_path), "Path to input single cell OH encoding csv does not exist."

        if self.verbose:
            print(f'Generating reference single cell input text files from {input_sc_tab_path}\nSaving results in directory: {directory}')

        subprocess.run(
            [
                path_to_C_file,
                input_sc_tab_path,
                input_labels_oh_path,
                directory,
                f"{n_topics}"
            ]
        )

        if self.verbose:
            print(f'Done generating reference GTM-formatted text data files to directory: {directory}')


    def generate_bulk_input_C_implementation(self, directory, path_to_C_file, input_bulk_tab_path):
        '''
        If bulk data is pre-formated in tab separated format, can use much quicker C implementation to generate bulkData.txt

        Args:
            directory : str
                - directory where GTM input text files will be saved
            path_to_C_file : str
                - path to C executable
            input_bulk_tab_path : str
                - path to bulk tab file

        Returns:
            None
        '''


        if self.verbose:
            print('Generating GTM format bulk text file...')

        subprocess.run(
            [
                path_to_C_file,
                input_bulk_tab_path,
                os.path.join(directory, BULKDATA_FILE),
            ]
        )

        if self.verbose:
            print(f'Done generating GTM format bulk text file at: {os.path.join(directory, BULKDATA_FILE)}')


    def generate_metagene_normalized(self, directory, inferenceMethod="JCVB0"):
        '''
        Generates metagene_normalized csv files, should be used only when n_topics != 1

        Args:
            directory : str
                - path to directory containing metagene files (deconvolved output)
            inferenceMethod : str
                - Don't change this
        
        Returns:
            None
        '''

        K = self.n_topics * len(self.celltypes)

        metagenes = os.listdir(directory)
        metagenes = [re.search(f"{BULKDATA_FILE[:-4]}_{REFERENCE_FILE[:-4]}_{inferenceMethod}_nmar_K{K}_iter\d+_metagene.csv", m) for m in metagenes]
        metagenes = sorted([m.group(0) for m in metagenes if m]) 

        if len(metagenes) == 0:
            warnings.warn("No metagene files were found, are you sure you deconvolved data to this directory?")
            return None

        for m in metagenes:
            metagene_df = pd.read_csv(os.path.join(directory, m), header=None)
            metagene_normalized = pd.DataFrame()

            for i in range(len(self.celltypes)):
                metagene_normalized[i] = metagene_df.iloc[:, i * self.n_topics : (i + 1) * self.n_topics].sum(axis=1)
        
            if self.verbose:
                print(f'Converting {m} to normalized metagene')

            metagene_normalized.to_csv(
                os.path.join(directory, m.replace('metagene', 'metagene_normalized')), 
                header=False, index=False
            )


    def gather_results(self, directory, output_path='', inferenceMethod="JCVB0"):
        '''
        Gathers metagene (if n_topics = 1) and metagene_normalized (n_topics > 1) into single results file, and labels celltype/batch information.

        Args:
            directory : str
                - directory containing metagene and metagene_normalized files (wherever deconvolved output was saved to)
            output_path : str [optional]
                - if you want to save to a different file name or directory
                - by default will save to directory with name gatheredResults.csv
            inferenceMethod : str
                - Don't change this
        
        Returns:
            metagene_df : pd.DataFrame
                - DataFrame containing gathered information if you want to use it immediately
        '''


        K = self.n_topics * len(self.celltypes)
        metagenes = os.listdir(directory)
        if self.n_topics == 1:
            metagenes = [re.search(f"{BULKDATA_FILE[:-4]}_{REFERENCE_FILE[:-4]}_{inferenceMethod}_nmar_K{K}_iter\d+_metagene.csv", m) for m in metagenes]
        else:
            metagenes = [re.search(f"{BULKDATA_FILE[:-4]}_{REFERENCE_FILE[:-4]}_{inferenceMethod}_nmar_K{K}_iter\d+_metagene_normalized.csv", m) for m in metagenes]
           
        metagenes = sorted([m.group(0) for m in metagenes if m]) 

        if len(metagenes) == 0:
            warnings.warn("No deconvolved results were found in this directory, check if there are any metagene/metagene_normalized files.")
            return None

        iter_list = []
        metagene_list = []
        for m in metagenes:
            metagene_list.append(pd.read_csv(os.path.join(directory, m), header=None))
            iter_list.append(re.sub("\D", '', re.search("iter\d+", m).group(0)))
        
        metagene_df = pd.concat(metagene_list)
        metagene_df.columns = self.celltypes

        # only save iter info if run for intermediates are saved
        if len(iter_list) != 1:
            l = [[i] * len(self.bulk_samples) for i in iter_list]
            l = [int(item) for sublist in l for item in sublist]
            metagene_df['iter'] = l

        metagene_df.index = self.bulk_samples * len(metagenes)


        
        # output_path is optional, default behavior is to dump to GATHERED_RESULTS_FILE in same directory as GTM output
        if output_path:
            if self.verbose:
                print(f'Compiling deconvolution results to {output_path}')
            metagene_df.to_csv(output_path)
        else:
            if self.verbose:
                print(f'Compiling deconvolution results to {os.path.join(directory, GATHERED_RESULTS_FILE)}')
            metagene_df.to_csv(os.path.join(directory, GATHERED_RESULTS_FILE))

        return metagene_df


    def pipeline(self, 
        bulk_data, 
        reference_data, 
        directory,
        ref_adata_var_gene_column=None,
        ref_adata_obs_celltype='Celltype', 
        ref_adata_layers_name=None, 
        bulk_adata_var_gene_column=None,
        bulk_adata_obs_sample_name=None,
        bulk_adata_layers_name=None, 
        n_topics=None, 
        engine_path=None, 
        n_iter=5,
        maxcores=10,
        output_intermediates=False,
        ):

        if self.verbose:
            print('Running GTM Deconvolution Pipeline')
            print(f'Writing results to {directory}')
            print('**********************************\n')

        if isinstance(bulk_data, pd.DataFrame):
            bulk_df, bulk_ad = bulk_data, None
        elif isinstance(bulk_data, anndata.AnnData):
            bulk_df, bulk_ad = None, bulk_data

        if isinstance(reference_data, pd.DataFrame):
            ref_df, ref_ad = reference_data, None
        elif isinstance(reference_data, anndata.AnnData):
            ref_df, ref_ad = None, reference_data


        if n_topics:
            self.n_topics = n_topics
        if engine_path:
            self.engine_path = engine_path
        
        self.generate_reference_input(
            directory=directory,
            reference_DataFrame=ref_df,
            reference_AnnData=ref_ad,
            adata_var_gene_column=ref_adata_var_gene_column,
            adata_obs_celltype=ref_adata_obs_celltype,
            adata_layers_name=ref_adata_layers_name,
        )
        
        self.generate_bulk_input(
            directory=directory,
            bulk_DataFrame=bulk_df,
            bulk_AnnData=bulk_ad,
            adata_var_gene_column=bulk_adata_var_gene_column,
            adata_obs_sample_name=bulk_adata_obs_sample_name,
            adata_layers_name=bulk_adata_layers_name,
        )

        self.save_parameters(directory)

        self.run_training(
            directory=directory,
            n_iter=n_iter,
            maxcores=maxcores,
            output_intermediates=output_intermediates,
        )

        self.run_deconvolution(
            directory=directory,
            maxcores=maxcores,
        )

        if self.n_topics != 1:
            self.generate_metagene_normalized(directory=directory)

        self.gather_results(directory=directory)


        if self.verbose:
            print('Completed running Pipeline.')



if __name__ == "__main__":

    pass


