#!python3

'''
Example Leave-One-Out Cross Validation Script using paired single cell reference and real bulk data.
'''

import pandas as pd
import os
import sys 

sys.path.append( os.path.dirname( os.path.dirname (os.path.abspath(__file__) ) ) )
from PythonWrapper.GTM_decon import GTM_decon





if __name__ == "__main__":
    # check if intended example data exists
    assert os.path.exists("../data/bulk_data.csv"), "Please first go through the tutorial Jupyter Notebook or extract ../data/tutorial_data.tar.gz"
    assert os.path.exists("../data/reference_data.csv"), "Please first go through the tutorial Jupyter Notebook or extract ../data/tutorial_data.tar.gz"


    bulk = pd.read_csv("../data/bulk_data.csv", index_col=0)
    reference = pd.read_csv("../data/reference_data.csv")
    batches = sorted(set(reference['Batch']))

    # if necessary, make results directory
    if not os.path.exists('tutorial_results'):
        os.mkdir('tutorial_results')

    # Run LOOCV in separate results sub-directories
    for batch in batches:
        if not os.path.exists(os.path.join('tutorial_results', batch)):
            os.mkdir(os.path.join('tutorial_results', batch))
        directory = os.path.join('tutorial_results', batch)

        reference_df = reference[reference['Batch'] != batch].reset_index(drop=True).drop(columns=['Batch'])
        bulk_df = bulk[[batch]]


        GTM = GTM_decon(
            engine_path = "/home/mcb/users/slaksh1/projects/revision_gb/gtm-decon-phinorm/gtm-decon-plus-noupd-ab-phinorm", #TODO
            experiment_name = f'vignette-leave-out-{batch}',
            n_topics = 5,
            verbose = False,
        )

        GTM.pipeline(
            directory = directory,
            reference_data = reference_df,
            bulk_data = bulk_df,
            n_iter=4
        )
    


