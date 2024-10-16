import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from beta_estimation_functions import *


def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

def count_total_elements(my_list):
    total_elements = 0

    for item in my_list:
        if isinstance(item, list):
            total_elements += count_total_elements(item)
        else:
            total_elements += 1

    return total_elements

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


def get_beta_all_func(beta_function_fid):

  if beta_function_fid == "@":
      def beta_function(P_mat, beta_arr):
          assert len(P_mat[0]) == len(beta_arr)
          P_mat = np.nan_to_num(np.array(P_mat), nan=1e-30)
          beta_arr = np.nan_to_num(np.array(beta_arr), nan=1e-30)
          P_mat[P_mat == 0] = np.float64(1e-30)
          beta_arr[beta_arr == 0] = np.float64(1e-30)
          return P_mat @  beta_arr.T

  elif beta_function_fid == "A":
      def beta_function(P_mat, beta_arr):
          assert len(P_mat[0]) == len(beta_arr)
          P_mat = np.nan_to_num(np.array(P_mat), nan=1e-30)
          beta_arr = np.nan_to_num(np.array(beta_arr), nan=1e-30)
          P_mat[P_mat == 0] = np.float64(1e-30)
          beta_arr[beta_arr == 0] = np.float64(1e-30)
          return np.log10(P_mat) @ np.log10(beta_arr).T

  elif beta_function_fid == "B":
    def beta_function(P_mat, beta_arr):
        P_mat = np.nan_to_num(np.array(P_mat, dtype=np.float64), nan=1e-30)
        beta_arr = np.nan_to_num(np.array(beta_arr, dtype=np.float64), nan=1e-30)
        
        P_mat[P_mat == 0] = 1e-30
        beta_arr[beta_arr == 0] = 1e-30
        
        return np.nan_to_num(P_mat) @ np.log10(beta_arr).T

  return beta_function

def tf_probability(Kd_tf, N_tf, genome_len):
  if tf_probabiltiy_fid == "@":
    return N_tf / (genome_len * (N_tf + Kd_tf))
  elif tf_probabiltiy_fid == "A":
    return N_tf / (N_tf + Kd_tf)

def get_coefficient_mat(current_protein_amnts, gene_info_dict, gene_key, tf_key):
  coefficient_mat = np.zeros((1, len(tf_key)))
  for tf, info in gene_info_dict["regulators"].items():
      if tf ==  "polymerase":
         continue
      Kd_tf_tg = score_to_K(info["delta G"], temperature)
      try:
        N_tf = current_protein_amnts[gene_key.index(tf)]
        coefficient_mat[0, tf_key.index(tf)] = tf_probability(Kd_tf_tg, N_tf, genome_length)
      except:
        continue
  return coefficient_mat

def get_R_mRNA_txn(N_rnap, R_max_txn, beta_all, Kd_rnap, genome_length, transcriptionRate_betaFromContext_fid):
  if transcriptionRate_betaFromContext_fid == "@":
      P_rnap = N_rnap / (genome_length * (N_rnap + Kd_rnap))
      return beta_all * (P_rnap * R_max_txn)
  elif transcriptionRate_betaFromContext_fid == "A":
      P_rnap = N_rnap / (N_rnap + Kd_rnap)
      return beta_all * (P_rnap * R_max_txn)
  elif transcriptionRate_betaFromContext_fid == "B":
      return  beta_all * (N_rnap * R_max_txn / (Kd_rnap + 1))

  # return np 1d arr of all txn rates


def mRNA_change(gene_mRNA_amnt, protein_amnts, gene_info_dict, N_rnap, R_max_txn):
  return get_R_mRNA_txn(
    N_rnap, R_max_txn, protein_amnts, gene_info_dict,
  ) - RNA_decay_rate(
     gene_mRNA_amnt, protein_amnts, gene_info_dict, gene_key, temperature, growthRate_fid, mRNA_decay_rate_farr
  ) * gene_mRNA_amnt


def protein_change(gene, gene_prot_amnt, protein_amnts, gene_info_dict, N_ribo, gene_mRNA_amnt, R_max_trans):
  return translation_rate(
     gene, gene_mRNA_amnt, protein_amnts, N_ribo, gene_info_dict, R_max_trans, Kd_ribo_mrna, temperature, translationRate_fid
    ) * gene_mRNA_amnt - protein_decay_rate(
       gene, protein_amnts, gene_info_dict, gene_key, temperature, growthRate_fid, protein_decay_rate_farr
    ) * gene_prot_amnt

def test_config_accuracy(testing_df, tf_key, beta_values, network_loc, network_key, feature_id):

  define_const_and_fids(feature_id)
  global ntf_key
  ntf_key = tf_key
  global gene_key
  global beta_arr
  beta_arr = beta_values

  gene_key = list(testing_df.index)
  time_points_key = list(testing_df.columns)
  actual_mRNAs_matrix = testing_df.to_numpy().T # TIME POINTS X GENE NAMES
  predicted_mRNAs_matrix = np.array(np.zeros((1, len(gene_key)))) # TIME POINTS X GENE NAMES: the first time point cannot be estimated so zeros

  global beta_func
  beta_func = get_beta_all_func()
  protein_t0_predict = actual_mRNAs_matrix[0, :]

  N_rnap = 2200
  N_ribo = 3000

  # for each indvidual time pair
  for t0 in range(len(time_points_key) - 1):
    mRNA_t0 = actual_mRNAs_matrix[t0, :]
    dt = int(time_points_key[t0]) - int(time_points_key[t0+1])

    mRNA_t1_predict = []
    protein_t1_predict = []

    for i in range(len(gene_key)):
      gene = gene_key[i]
      mRNA_t0_gene = mRNA_t0[i]
      prot_t1_gene = protein_t0_predict[i]

      gene_info_dict = search_network_json(network_loc, network_key, gene)
      if gene_info_dict == None:
         protein_t1_predict.append(None)
         mRNA_t1_predict.append(np.float64("nan"))
         continue
      R_max_txn, R_max_trans = max_rates(len_taken_by_rnap, elongation_rate, len_taken_by_ribo, peptide_rate, gene_info_dict, transcriptionRate_betaFromContext_fid)

      # get mrna amounts
      dmRNAdt = mRNA_change(mRNA_t0_gene, protein_t0_predict, gene_info_dict, N_rnap, R_max_txn)
      gene_mRNA_t1_predict = mRNA_t0_gene + dmRNAdt * dt
      gene_mRNA_t1_predict = gene_mRNA_t1_predict * 6.022140857e23 / 1e15 # concentration to amounts
      mRNA_t1_predict.append(gene_mRNA_t1_predict)

      # determine the protein amnts for next round
      dPdt = protein_change(gene, prot_t1_gene, protein_t0_predict, gene_info_dict, N_ribo, mRNA_t0_gene, R_max_trans)
      protein_t1_predict_gene = prot_t1_gene + dPdt * dt
      protein_t1_predict.append(protein_t1_predict_gene)

    mRNA_t1_predict = np.array(mRNA_t1_predict, "object")
    predicted_mRNAs_matrix = np.vstack([predicted_mRNAs_matrix, mRNA_t1_predict])
    protein_t0_predict = protein_t1_predict

    N_rnap = RNAP_amount(N_rnap, protein_t0_predict, RNAPAmount_fid)
    N_ribo = ribo_amount(N_ribo, protein_t0_predict, RiboAmount_fid)

  return actual_mRNAs_matrix, predicted_mRNAs_matrix, gene_key, time_points_key

# get accuracy: compare predict to actual
def accuracy_representation(config_id, data_actual_matrix, data_predicted_matrix, gene_key, time_points_key, all_bad_indecies_tup, prokode_dir, just_overall):
  print("AA")
  flattened_all_bad_tup = [x for xs in all_bad_indecies_tup for x in xs]
  flattened_actual = data_actual_matrix.flatten()[flattened_all_bad_tup]
  flattened_predict = data_predicted_matrix.flatten()[flattened_all_bad_tup]
  print(flattened_predict)

  if just_overall:

    print("dk")
    error_raw = np.nan_to_num(flattened_actual) - np.nan_to_num(flattened_predict)
    print('a', np.isnan(np.nan_to_num(flattened_actual)).any())
    print('b', np.isnan(np.nan_to_num(flattened_predict)).any())

    print("raw", error_raw)
    print(np.isnan(error_raw).any())
    mse_error = np.sum(np.nan_to_num(np.square(error_raw))) / len(flattened_actual)
    return mse_error

  # create normal probability plot
  plt.scatter(data_actual_matrix.flatten(), data_predicted_matrix.flatten())
  plt.xlabel('Actual gene expression (mRNA in mol / L)')
  plt.ylabel('Predicted gene expression (mRNA in mol / L)')
  plt.axline((0,0), (1,1))
  plt.title('Actual vs Predicted Petal Width')
  plt.savefig(f'{prokode_dir}/beta_training/config_results/{config_id}/normal_prob.pdf')

  # maybe area under curve or smth
  test_df = pd.DataFrame(np.concat((flattened_actual, flattened_predict), axis=1), index = time_points_key, columns=["actual", "predict"])
  test_df.to_csv(f"{prokode_dir}/beta_training/config_results/{config_id}/raw_predic_act.csv")

  return "Success!"
