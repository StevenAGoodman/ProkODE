import numpy as np
import pandas as pd
import scipy as sp
import os
from run import *

prokode_dir = '/workspaces/PROKODE-DOCKER'

def get_beta_all_static():

def get_beta_all_timeseries():
    
    beta_all = 
    return beta_all

def get_betas(network_json, sample_gene_dict tf_key):
    coefficient_matrix = []
    
    # create coefficient array by adding the P_tg things then create beta_all matrix 
    for tg in sample_gene_dict:
        tg_info = network_json[tg]
        # define tfs
        tfs_info = [tf for tf in tg_info["regulators"]]

        coefficient_arr = [0 for i in len(tf_key)]
        
        # find P_m-n s
        for tf in tfs_info:
            N_tf = sample_gene_dict[tf]
            P = N_tf / (N_tf + tf["kd"])
            coefficient_arr[tf_key.index(tf)] = P
            
        coefficent_matrix.append(coefficient_arr)

        # beta_all_matrix
        beta_all = get_beta_all()
        beta_all_matrix.append(beta_all)
    
    coefficient_matrix = np.array(coefficient_matrix)
    beta_all_matrix = np.array(beta_all_matrix, columnshape)
    
    _, beta_arr = sp.linalg.sparse.lsqr(coefficient_matrix, beta_all_matrix)
    return beta_arr

def get_betas(sample_loc, genome_loc, annotation_loc):
    # read sample to df
    sample_df = pd.read_csv(sample_loc, names=['gene', 'transcription rate'])

    net_json_loc = create_network_json_main(prokode_dir,)
    network_dict = json.load(open(net_json_loc, 'r'))

    for _,row in sample_df.iter_rows:
        gene = row[0]
        transcription_rate = row[1]

        get_betas(transcription_rate, N_tf, Kd_tf, N_p, Kd_p, [all tf jazz]) # from maths.py get beta value from context & trans rate

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
