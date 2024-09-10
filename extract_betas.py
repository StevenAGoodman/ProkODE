import numpy as np
import pandas as pd
import os


def get_betas():
    # read sample to df
    sample_df = pd.read_csv(sample_file, names=['gene', 'transcription rate'])

    fit_beta

    return tf_betas

def main(data_folder_loc, results_loc):
    # import all expiremental data files into folder
    data_folder = os.fsencode(data_folder_loc)
    
    results_csvoutput = ''

    # loop over all expiremental sample files in folder
    for sample_file in os.listdir(data_folder):
        sample_loc = os.fsdecode(sample_file)

        # get str in csv formatting of tf - beta relations for sample
        sample_tf_betas = get_betas()

        results_csvoutput += sample_tf_betas

    # write results to training csv (tf, beta)
    open(results_loc, 'w').write(results_csvoutput)

main('/data/samples')