import pandas as pd
import numpy as np
import json
import statistics
from src.run import create_network_json_main
import cvxpy

# be sure to run run.py before running this file
network_loc = "src/network.json"
data_loc = "./GSE10159_results.csv"

# global jazz
num_data = 10
sensor_normal_dist = 10
basal_rate = 3
decay_rate = np.log(2)/300
Np = 6000
Kd_p = 0.1
Nns = 4600000

# def function_to_test(data_t0, protein_data_t0, dt):
     # rate_of_protein_creation = some_function_dependent_on_num_ribo()
     # overall_protein_half_life = some_constant_func_independent_of_other_stuff()

     # overall_rate_of_protein_change = rate_of_protein_creation * data_t0  - overall_protein_half_life * protein_data_t0

     # rate_of_mRNA_creation =
     # overall_mRNA_half_life = some_constant_func_independent_of_other_stuff()

     # overall_rate_of_mRNA_change = rate_of_mRNA_creation - overall_mRNA_half_life * data_t0

     # protein_data_t1 = protein_data_t0 + overall_rate_of_protein_change * dt
     # data_t1 = data_t0 + overall_rate_of_mRNA_change * dt
     # return data_t1

def protein_amnts_from_mRNA_amnts(mRNA_amnts):
     protein_amnts = mRNA_amnts
     return protein_amnts

from scipy.optimize import curve_fit
from src.maths import *

def fit_function_for_timepoint(data_t0, data_t1, dt, gene_key):
     # MOST CRITICAL PART TO DETERMINE!!!
     protein_data_t0 = protein_amnts_from_mRNA_amnts(data_t0)

     for i in range(len(gene_key)):
          gene = gene_key[i]
          gene_mRNA_t0 = data_t0[i]

          overall_mRNA_change_rate = (data_t1[i] - data_t0[i]) / dt
          beta_all_arr, coefficient_matrix = creation_rate_from_overall_mRNA_stats(overall_mRNA_change_rate, gene_mRNA_t0, protein_data_t0, gene_key)

          def beta_model_func(x_arr, beta_arr):
               y = to_beta_all(x_arr, beta_arr)
               return y

          beta_arr = curve_fit(beta_model_func, xdata = coefficient_matrix, ydata = beta_all_arr)[0]

     return(beta_arr)



# input: pair of time points

data_df = pd.read_csv(data_loc, index_col=0)
gene_key = list(data_df.index) # may have to subtracted header

# lil preprocessing


for i in range(len((data_df.index) - 1)):
     data_t0 = list(data_df.iloc[:,i])
     data_t1 =  list(data_df.iloc[:,i + 1])
     dt = int(data_df.columns[i + 1]) - int(data_df.columns[i])

     beta_arr = fit_function_for_timepoint(data_t0, data_t1, dt, gene_key)



