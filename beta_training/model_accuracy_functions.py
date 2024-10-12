import numpy as np
import pandas as pd

def test_config_accuracy(testing_df, tf_key, beta_values):

  gene_key = list(testing_df.columns())
  time_points_key = list(testing_df.indecies())

  actual_mRNAs_matrix = testing_df.to_numpy() # TIME POINTS X GENE NAMES

  predicted_mRNAs_matrix = np.array(np.zeros((1, len(gene_key)))) # TIME POINTS X GENE NAMES: the first time point cannot be estimated so zeros

  # for each indvidual time pair
  for t0 in range(len(testing_df) - 1):
    mRNA_t0 = testing_df[:, t0]
    dt = int(testing_df.columns[t0]) - int(testing_df.columns[t0+1])

    # get mrna amounts
    dmRNAdt = mRNA_change() #np arr
    mRNA_t1_predict = mRNA_t0 + dmRNAdt * dt

    predicted_mRNAs_matrix = np.concat((predicted_mRNAs_matrix, mRNA_t1_predict),axis=0)


    # determine the protein amnts for next round
    dPdt = protein_change()
    protein_t1_predict = protein_t0_predict + dPdt * dt
    protein_t0_predict = protein_t1_predict

  return actual_mRNAs_matrix, predicted_mRNAs_matrix, gene_key, time_points_key

# get accuracy: compare predict to actual
def measure_accuray(data_actual_matrix, data_predicted_matrix, gene_key, time_points_key):
  # create predict vs actual graph
  # measure maean squared accuracy
  # maybe area under curve or smth
  return some_measure of accuracy
