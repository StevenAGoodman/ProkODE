import sys
# inputs
feature_string = sys.argv[1] # "@@@@@@@@10001000" #
data_file = sys.argv[2] # 
prokode_dir = sys.argv[3] # "/workspaces/ProkODE" # 

import pandas as pd
import numpy as np
import json
import statistics
import cvxpy
import re
import os
import multiprocessing as mp
 
sys.path.append(".")
from scipy.optimize import curve_fit
from maths import *

# be sure to run run.py before this file
network_loc = prokode_dir + "/src/network.json"
data_dir = prokode_dir + "/beta_training/GEO_expression_data" 
# data_file_list = os.listdir(data_dir)

# constants that preferably would be set by the user
genome_length = 4500000
cell_volume = 1e-15 # L
temperature = 298 # kelvin
elongation_rate = 60 # nt/s
peptide_rate = 20 # aa/s
Kd_ribo_mrna = 0.1 # WHAT IS THE BINDING AFFINITY OF RIBO TO DALGARNO SEQ???? maybe use molecular docking
len_taken_by_rnap = 30 # nt
len_taken_by_ribo = 30 # aa

# unpack model feature string
feature_ids = list(feature_string)
transcriptionRate_betaFromContext_fid = feature_ids[0]
translationRate_fid = feature_ids[1]
tf_probabiltiy_fid = feature_ids[2]
beta_function_fid = feature_ids[3]
RNAPAmount_fid = feature_ids[4]
RiboAmount_fid = feature_ids[5]
sigma_competition_fid = feature_ids[6]
growthRate_fid = feature_ids[7]
mRNA_decay_rate_farr = feature_ids[8:12]
protein_decay_rate_farr = feature_ids[12:16]

# global definable functions
network_key = create_network_key(network_loc)
func_to_fit = beta_curveFit_func(beta_function_fid)

def function_for_timepoint(mRNAs_t0, N_rnap, N_ribo, protein_data_t0, cell_volume, gene_key, tf_key):
     protein_data_t1 = [None] * len(gene_key)
     coefficient_submat = []
     beta_all_subarr = []

     for i in range(len(gene_key)):
          gene = gene_key[i]
          gene_mRNA_t0 = data_t0[i]
          gene_info_dict = search_network_json(network_loc, network_key, gene)

          if gene_info_dict == None:
               print(f"\t\t\t\tFailed: {gene}\tgene from data not found in genome")
               continue
          
          R_max_txn, R_max_trans = max_rates(len_taken_by_rnap, elongation_rate, len_taken_by_ribo, peptide_rate, gene_info_dict, transcriptionRate_betaFromContext_fid)
          overall_mRNA_change_rate = (data_t1[i] - data_t0[i]) / dt

          coefficient_arr, beta_all = beta_from_overall_mRNA(
               overall_mRNA_change_rate, 
               gene_mRNA_t0,
               protein_data_t0,
               N_rnap, 
               gene_info_dict,
               gene_key,
               tf_key,
               R_max_txn,
               temperature, genome_length,
               growthRate_fid, tf_probabiltiy_fid, sigma_competition_fid, transcriptionRate_betaFromContext_fid, mRNA_decay_rate_farr
          )

          coefficient_submat.append(coefficient_arr)
          beta_all_subarr.append(beta_all)

          # update protein amnts
          overall_protein_change_rate = translation_rate(
               gene, gene_mRNA_t0, protein_data_t0, N_ribo, gene_info_dict, R_max_trans, Kd_ribo_mrna, temperature, translationRate_fid
          ) - protein_decay_rate(
               gene, protein_data_t0, gene_info_dict, gene_key, temperature, growthRate_fid, protein_decay_rate_farr
          ) * protein_data_t0[i]

          protein_data_t1[i] = protein_data_t0[i] + overall_protein_change_rate * dt

     # update protein amounts, RNAP, ribo, and cell vol
     next_N_rnap = RNAP_amount(N_rnap, protein_data_t0, RNAPAmount_fid)
     next_N_ribo = ribo_amount(N_ribo, protein_data_t0, RiboAmount_fid)
     cell_volume = cell_volume + get_grow_rate(protein_data_t0, growthRate_fid) * dt


     return coefficient_submat, beta_all_subarr, N_rnap, N_ribo, protein_data_t1, cell_volume

def fit_function(coefficient_matrix, beta_all_arr, feature_id):
     # clean matrices
     beta_all_arr = np.nan_to_num(np.array(beta_all_arr)).tolist()

     empty_columns = []
     for i in range(len(list(zip(*coefficient_matrix)))):
          column = list(zip(*coefficient_matrix))[i]
          if column == [0] * len(coefficient_matrix):
               empty_columns.append(i)
          else:
               continue

     n_features = len(coefficient_matrix[0])
     p = [1] * n_features
     try:
          beta_arr, covar_uncertainty = curve_fit(func_to_fit, coefficient_matrix, beta_all_arr, p0 = p)

          for index in empty_columns:
               beta_arr[index] = None

          return beta_arr.tolist(), covar_uncertainty.tolist()
     except:
          return "ERROR IN FITTING", "NA"

for data_file in [data_file]:
     # get data
     data_df = pd.read_csv(f"{data_file}", index_col="Unnamed: 0").iloc[:, :3].multiply(1 / 6.02e23)
     print(data_df)
     # normalize as concentations
     # data_df = data_df.multiply(1 / 6.02e23)

     # define keys
     gene_key = list(data_df.index) # may have to subtracted header
     tf_key = [ tf for val in json.load(open(network_loc, 'r')).values() for tf in val["regulators"].keys() ]
     tf_key = list(set(tf_key)) # filter repeats

     # filter out groups
     groups = [ [] for i in range(1000) ]
     for i in range(len(data_df.axes[1])):
          colname = data_df.columns[i]
          group_n = re.search("\\| (\\d+)", colname)
          if group_n == None:
               group_n = "0"
          else:
               data_df = data_df.rename(columns = {colname: colname[:colname.find(" | ")]})
               group_n = group_n.group(1)

          groups[int(group_n)].append(i)
          
     groups = [x for x in groups if x != []] # filter empties

     beta_all_arr, coefficient_matrix = ([], [])

     # groups are isolated consecutive data
     for i in range(len(groups)):
          print(f"\tgroup {i}")

          group_data_indecies = sorted(groups[i])

          protein_data_t0 = [ float(i) / cell_volume for i in list(data_df.iloc[:,0]) ]
          N_rnap = 5000
          N_ribo = 15000

          for l in range(len(group_data_indecies[:-1])):
               print(f"\t\tdata pair {l}")

               n = group_data_indecies[l]  # data index at t0
               m = group_data_indecies[l+1]  # data index at t1

               # pair of time points
               data_t0 = [ float(i) / cell_volume for i in list(data_df.iloc[:,n]) ]  # data  at t0
               data_t1 = [ float(i) / cell_volume for i in list(data_df.iloc[:,m]) ]  # data  at t1
               dt = float(data_df.columns[m]) - float(data_df.columns[n])

               coefficient_submat, beta_all_subarr, N_rnap, N_ribo, protein_data_t0, cell_volume = function_for_timepoint(data_t0, N_rnap, N_ribo, protein_data_t0, cell_volume, gene_key, tf_key)
               
               coefficient_matrix.extend(coefficient_submat)
               beta_all_arr.extend(beta_all_subarr)

               # print(f"\t\t\tbeta all: {beta_all_subarr}")

# overall function fitting for all data points in an organisms
beta_arr, covar = fit_function(coefficient_matrix, beta_all_arr, beta_function_fid)

# export to csv
results_dict = {feature_string:{
     "beta values": dict(zip(tf_key, beta_arr)),
     "covariance": str(covar)
}}
results_dictstr = json.dumps(results_dict)
open(prokode_dir +"/beta_training/tf_betas.json", 'a').write(results_dictstr + "\n")