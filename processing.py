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

def get_data_unit(i):

def fit_function_for_timepoint(data_item, gene_key):
     data_t0 = data_item[:,0]
     data_t1 = data_item[:,1]
     

     return(beta_arr)
    
    
# input: pair of time points

data_df = pd.read_csv(data_loc)
gene_key = data_df[0] # may have to subtracted header

for i in range(names(data_df)):
    data_unit = get_data_unit(i)

    
