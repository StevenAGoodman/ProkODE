import numpy as np
import pandas as pd
import os
from run import *

def network_hub():

def get_betas(sample_loc):
    # read sample to df
    sample_df = pd.read_csv(sample_loc, names=['gene', 'transcription rate'])

    tf_ref_arr, = network_hub()
    
    for _,row in sample.iter_rows:
        gene = row[0]
        transcription_rate = row[1]
        
        get_gene_info()

        beta_from_eqs(transcription_rate, N_tf, Kd_tf, N_p, Kd_p, [all tf jazz]) # from maths.py get beta value from context & trans rate

    return tf_betas

def preprocessing_main(data_folder_loc, results_loc):
    # import all expiremental data files into folder
    data_folder = os.fsencode(data_folder_loc)
    
    results_csvoutput = ''

    # loop over all expiremental sample files in folder
    for sample_file in os.listdir(data_folder):
        sample_loc = os.fsdecode(sample_file)

        # get str in csv formatting of tf - beta relations for sample
        sample_tf_betas = get_betas(sample_loc)

        results_csvoutput += sample_tf_betas

    # write results to training csv (tf, beta)
    open(results_loc, 'w').write(results_csvoutput)

preprocessing_main('/data/samples')
