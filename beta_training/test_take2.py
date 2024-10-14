NUM_CORES = 4
PYTHON_EXEC_PATH = "python"
TF_BETA_FILEPATH = "c:/Users/cryst/LOFScreening/archive/PROKODE/ProkODE/partially_run_results.json"
DATA_FILE_PATH = "c:/Users/cryst/LOFScreening/archive/PROKODE/ProkODE/beta_training/GEO_expression_data/combined_gene_expression.csv"
PROKODE_BASEDIR =  "c:/Users/cryst/LOFScreening/archive/PROKODE/ProkODE" # !! NO end '/' !!

# libs
import multiprocessing as mp
import time
import os
import pandas as pd
import numpy as np
from beta_estimation_functions import *
from model_accuracy_functions import *

def define_const_and_fids(feature_string):
  # constants that preferably would be set by the user
  global genome_length
  genome_length = 4500000
  global cell_volume
  cell_volume = 1e-15 # L
  global temperature
  temperature = 298 # kelvin
  global elongation_rate
  elongation_rate = 60 # nt/s
  global peptide_rate
  peptide_rate = 20 # aa/s
  global Kd_ribo_mrna
  Kd_ribo_mrna = 0.1 # WHAT IS THE BINDING AFFINITY OF RIBO TO DALGARNO SEQ???? maybe use molecular docking
  global len_taken_by_rnap
  len_taken_by_rnap = 30 # nt
  global len_taken_by_ribo
  len_taken_by_ribo = 30 # aa

  # unpack model feature string
  feature_ids = list(feature_string)
  global transcriptionRate_betaFromContext_fid
  transcriptionRate_betaFromContext_fid = feature_ids[0]
  global translationRate_fid
  translationRate_fid = feature_ids[1]
  global tf_probabiltiy_fid
  tf_probabiltiy_fid = feature_ids[2]
  global beta_function_fid
  beta_function_fid = feature_ids[3]
  global RNAPAmount_fid
  RNAPAmount_fid = feature_ids[4]
  global RiboAmount_fid
  RiboAmount_fid = feature_ids[5]
  global sigma_competition_fid
  sigma_competition_fid = feature_ids[6]
  global growthRate_fid
  growthRate_fid = feature_ids[7]
  global mRNA_decay_rate_farr
  mRNA_decay_rate_farr = feature_ids[8:12]
  global protein_decay_rate_farr
  protein_decay_rate_farr = feature_ids[12:16]

  global Kd_rnap
  Kd_rnap = 10.1


def getGroups(data_df):
    # filter out groups
    print("\t\t\t\tfiltering different groups")
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
    minim = min(min(groups))
    groups = [[a - minim for a in r] for r in groups]

    return data_df, groups

def init_stuffs(data_file, network_loc):
    print("reading datafile to dataframe")
    # get data
    data_df = pd.read_csv(f"{data_file}", index_col="Unnamed: 0").multiply(1 / 6.02e23)
    # normalize as concentations
    # data_df = data_df.multiply(1 / 6.02e23)
    # define keys
    gene_key = list(data_df.index) # may have to subtracted header
    tf_key = [ tf for val in json.load(open(network_loc, 'r')).values() for tf in val["regulators"].keys() ]
    tf_key = list(set(tf_key)) # filter repeats

    data_df = data_df.iloc[:, 40:66]
    data_df, groups = getGroups(data_df)

    group_drops = []
    for m in range(len(groups)):
        group = groups[m]
        if len(group) < 2:
            group_drops.append(m)
    groups = [i for j, i in enumerate(groups) if j not in group_drops]
    return data_df, gene_key, tf_key, groups

def get_coefficient_arr(gene_info_dict, protein_amnts):
    coefficient_arr = [0] * len(tf_key)
    for regulator, reg_details in gene_info_dict["regulators"].items():
        if regulator == "polymerase":
            continue
        try:
            N_tf = protein_amnts[gene_key.index(regulator)]
            Kd_tf_target = score_to_K(reg_details["delta G"], temperature)
            P_tf = tf_probabiltiy(N_tf, Kd_tf_target, genome_length, tf_probabiltiy_fid)
            if np.isnan(N_tf) or np.isnan(Kd_tf_target) or  np.isnan(P_tf):
                continue
        except:
            # print(f"\t\t\tFailed: {regulator}\tgenome tf is not found in data's genes")
            continue
        coefficient_arr[tf_key.index(regulator)] = P_tf

    return [coefficient_arr]


def get_txn_and_prots_and_rnadecay(data_t0, protein_data_t0, dt, N_ribo, N_rnap, beta_values_arr, cell_volume):
    protein_data_t1 = [None] * len(gene_key)

    R_mRNA_txn_arr = []
    mRNA_decay_arr = []
    
    for i in range(len(gene_key)):
        gene = gene_key[i]
        gene_mRNA_t0 = data_t0[i]
        gene_info_dict = search_network_json(network_loc, network_key, gene)
        if gene_info_dict == None:
            # print(f"\t\t\t\tFailed: {gene}\tgene from data not found in genome")
            mRNA_decay_arr.append(np.nan)
            R_mRNA_txn_arr.append(np.nan)

            continue

        R_max_txn, R_max_trans = max_rates(len_taken_by_rnap, elongation_rate, len_taken_by_ribo, peptide_rate, gene_info_dict, transcriptionRate_betaFromContext_fid)
        coefficient_arr = get_coefficient_arr(gene_info_dict, protein_data_t0)

        # update protein amnts
        overall_protein_change_rate = translation_rate(
            gene, gene_mRNA_t0, protein_data_t0, N_ribo, gene_info_dict, R_max_trans, Kd_ribo_mrna, temperature, translationRate_fid
        ) - protein_decay_rate(
            gene, protein_data_t0, gene_info_dict, gene_key, temperature, growthRate_fid, protein_decay_rate_farr
        ) * protein_data_t0[i]
        protein_data_t1[i] = protein_data_t0[i] + overall_protein_change_rate * dt

        # what ya need
        beta_all = beta_all_func(coefficient_arr, beta_values_arr)
        R_mRNA_txn = get_R_mRNA_txn(N_rnap, R_max_txn, beta_all, Kd_rnap, genome_length, transcriptionRate_betaFromContext_fid) #check

        R_mRNA_txn_arr.append(R_mRNA_txn.tolist()[0])

        mRNA_decay_rate = RNA_decay_rate(gene_mRNA_t0, protein_data_t0, gene_info_dict, gene_key, temperature, growthRate_fid, mRNA_decay_rate_farr)
        mRNA_decay_arr.append(mRNA_decay_rate)

    # update protein amounts, RNAP, ribo, and cell vol
    next_N_rnap = RNAP_amount(N_rnap, protein_data_t0, RNAPAmount_fid)
    next_N_ribo = ribo_amount(N_ribo, protein_data_t0, RiboAmount_fid)
    cell_volume = cell_volume + get_grow_rate(protein_data_t0, growthRate_fid) * dt

    return R_mRNA_txn_arr, protein_data_t1, mRNA_decay_arr, next_N_rnap, next_N_ribo, cell_volume

def new_mRNA_change(R_mRNA_txn_arr, mRNA_decay_arr, data_t0_actual, dt):
    bad_indecies:list = np.where(np.isnan(R_mRNA_txn_arr))[0].tolist()
    bad_data_indecies:list = np.where(np.isnan(data_t0_actual))[0].tolist()
    bad_data_decay:list = np.where(np.isnan(mRNA_decay_arr))[0].tolist()
    all_bads:list = bad_indecies  # Start with the bad indices from R_mRNA_txn_arr
    all_bads.extend(bad_data_indecies)  # Add bad indices from data_t0_actual
    all_bads.extend(bad_data_decay)  # Add bad indices from mRNA_decay_arr
    
    data_t0_actual = np.nan_to_num(data_t0_actual, nan=0.0, posinf=1e-10, neginf=1e-10)
    R_mRNA_txn_arr = np.nan_to_num(R_mRNA_txn_arr, nan=0.0, posinf=1e-10, neginf=1e-10)
    mRNA_decay_arr = np.nan_to_num(mRNA_decay_arr, nan=0.0, posinf=1e-10, neginf=1e-10)
    
    dRdt_arr = R_mRNA_txn_arr.T - mRNA_decay_arr.T @ np.diag(data_t0_actual)
    dRdt_arr[all_bads] = data_t0_actual[all_bads] * dt
    dRdt_arr = dRdt_arr.tolist()
    # print('b', np.count_nonzero(np.isnan(dRdt_arr)))

    return dRdt_arr, all_bads


def main(testing_df, tf_key, beta_values_arr, network_loc, network_key, feature_string, groups):
    # be sure to run run.py before this file
    global all_bad_indecies_tup
    all_bad_indecies_tup = []

    data_predict_matrix = [] # times x genes
    gene_key = list(testing_df.index)
    time_points_key = list(testing_df.columns)
    
    # constants that preferably would be set by the user
    define_const_and_fids(feature_string)
    global beta_all_func
    beta_all_func = get_beta_all_func(beta_function_fid)
    global cell_volume

    # groups are isolated consecutive data
    for i in range(len(groups)):
        # print(f"\tgroup {i}")
        group_data_indecies = sorted(groups[i])
        protein_data_t0 = [ float(i) / cell_volume for i in list(testing_df.iloc[:,0]) ]
        N_rnap = 5000
        N_ribo = 15000
        for l in range(len(group_data_indecies[:-1])):
            # print(f"\t\tdata pair {l}")
            n = group_data_indecies[l]  # data index at t0
            m = group_data_indecies[l+1]  # data index at t1

            # pair of time points
            data_t0_actual = [ float(i) / cell_volume for i in list(testing_df.iloc[:,n]) ]  # data  at t0
            dt = float(testing_df.columns[m]) - float(testing_df.columns[n])

            # print("\t\t\testimating equs for each geme")
            R_mRNA_txn_arr, protein_data_t1, mRNA_decay_arr, N_rnap, N_ribo, cell_volume = get_txn_and_prots_and_rnadecay(
                data_t0_actual, protein_data_t0, dt, N_ribo, N_rnap, beta_values_arr, cell_volume
            )
            # print(R_mRNA_txn_arr)

            dmRNAdt, all_bad_indecies = new_mRNA_change(np.array(R_mRNA_txn_arr), np.array(mRNA_decay_arr), np.array(data_t0_actual), dt)
            # print(dmRNAdt[64:84])
            data_t1_predict = np.nan_to_num(data_t0_actual) + [ o * dt for o in dmRNAdt]
            # print(data_t1_predict[43:63])

            data_predict_matrix.append(data_t1_predict)
            protein_data_t0 = protein_data_t1
            all_bad_indecies_tup.append(all_bad_indecies)
            # print(f"\t\t\tbeta all: {beta_all_subarr}")

    data_actual_matrix = testing_df.to_numpy().T # time points x genes
    data_predict_matrix = np.array(data_predict_matrix) * (1e-15)

    return data_actual_matrix, data_predict_matrix, gene_key, time_points_key

   

network_loc = PROKODE_BASEDIR + "/src/network.json"
testing_df, gene_key, tf_key, groups = init_stuffs(DATA_FILE_PATH, PROKODE_BASEDIR + "/src/network.json")
print(testing_df, groups)
max_lines = 2594
network_key = create_network_key(network_loc)


def gauage(line_no):
    print("line:", line_no)

    line = linecache.getline(TF_BETA_FILEPATH, line_no)
    line = json.loads(line)

    beta_values = list([*line.values()][0]["beta values"].values())
    if type(beta_values[0]) == str:
        open("./beta_training/config_results/test_completed.txt", "a").write(str(line_no) + "\n")
        return "NA"     

    config_id = [*line.keys()][0]
    tf_key = list([*line.values()][0]["beta values"].keys())

    actual_matrix, predicted_matrix, gene_key, time_points_key = main(testing_df, tf_key, beta_values, network_loc, network_key, config_id, groups)
    print(actual_matrix.shape, predicted_matrix.shape)
    # do accuracy measurements
    general_acc = accuracy_representation(config_id, actual_matrix, predicted_matrix, gene_key, time_points_key, all_bad_indecies_tup, PROKODE_BASEDIR, True)
    if general_acc == np.nan:
        general_acc = str(general_acc) + ": floating point likely to small"
        

    open("./beta_training/config_results/overall.csv", 'a').write(f"{config_id},{general_acc}\n")
    open("./beta_training/config_results/test_completed.txt", "a").write(str(line_no) + "\n")

    return "Success!"



if __name__ == '__main__':
    pool = mp.Pool(processes=NUM_CORES)
    results = pool.map(gauage, list(range(max_lines)))
    pool.close()
    pool.join()
