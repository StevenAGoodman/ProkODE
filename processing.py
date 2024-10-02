import pandas as pd
import numpy as np
import json
import statistics
from src.run import create_network_json_main
import cvxpy 

# be sure to run run.py before running this file
network_loc =
data_loc =

# global jazz
num_data = 10
sensor_normal_dist = 10
basal_rate = 3
decay_rate = np.log(2)/300
Np = 6000
Kd_p = 0.1
Nns = 4600000

def function_to_test(data_t0, protein_data_t0, dt):
     rate_of_protein_creation = some_function_dependent_on_num_ribo()
     overall_protein_half_life = some_constant_func_independent_of_other_stuff()
     
     protein_data_t1 = rate_of_protein_creation * data_t0  - overall_protein_half_life * protein_data_t0

     rate_of_mRNA_cration = 
     overall_mRNA_half_life = some_constant_func_independent_of_other_stuff()
     
     overall_rate_of_mRNA_change = rate_of_mRNA_creation - overall_mRNA_half_life * data_t0

     data_t1 = data_t0 + overall_rate_of_mRNA_change * dt
     return data_t1

def fit_function_for_timepoint(data_t0, data_t1, gene_key):
     for i in range(len(gene_key)):
          gene = gene_key[i]

          rate_of_change
     
     return(beta_arr)
    
# input: pair of time points

data_df = pd.read_csv(data_loc)
gene_key = data_df[0] # may have to subtracted header

for i in range(len(names(data_df[1:] - 1))):
     data_t0 = data_df[i]
     data_t1 = data_df[i + 1])

     beta_arr = fit_function_for_timepoint(data_t0, data_t1, gene_key)
     

    
