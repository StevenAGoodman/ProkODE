import pandas as pd
import numpy as np
from scipy.sparse.linalg import lsqr
import json
import statistics
from run import create_network_json_main

prokode_dir = '/workspaces/PROKODE-DOCKER'
data_file = 'GSE90743_E14R025_raw_counts.txt'
genome_loc = '/workspaces/PROKODE-DOCKER/src/inputs/genome.fasta'
annotation_loc = '/workspaces/PROKODE-DOCKER/src/inputs/annotation.csv'

sensor_normal_dist = 10
basal_rate = 3
decay_rate = np.log(2)/300
Np = 6000
Nns = 4600000

def compare(beta_collection, tf_key):
    results = 'tf,mean,median,range,standard deviation\n'

    for tf_beta_index in range(len(tf_key)):
        tf = tf_key[tf_beta_index]
        tf_beta_vector = beta_collection[:, tf_beta_index]

        mean = statistics.fmean(tf_beta_vector)
        median = median(tf_beta_vector)
        range = [min(tf_beta_vector), max(tf_beta_vector)]
        std_dev = statistics.stdev(tf_beta_vector)

        results += f'{tf},{mean},{median},{range},{std_dev}\n'

    open('results.csv', 'w').write(results)

def update(P_predict, H, R, z_k, x_predict):
    # kalman gain
    K_k = P_predict @ H.transpose() @ np.linalg.inv(H @ P_predict @ H.transpose() + R)
    
    # update
    x_est = x_predict + K_k @ (z_k - H @ x_predict)
    P_k = P_predict - K_k @ H @ P_predict

    return x_est, P_k

def predict(x_prevk, P_prevk, A, Q):
    x_predict = A @ x_prevk
    P_predict = (A @ P_prevk @ A.transpose()) + Q
    return x_predict, P_predict

def kalman_filtering(data):
    gene_key = data.iloc[0, :].values.flatten().tolist()

    # define t
    t = [5,10,25,45,75,120,210,330,1500,1560,1680]

    # init kalman inputs
    x = x_0 = np.array([[0],
                  [0]])
    P = P_0 = np.array([[1000, 0],
                           [0, 1000]])
    F = np.array([[1, 1],
                  [0, 1]])
    Q = 0.01 # world error
    H = np.array([[1,0]])
    R = sensor_normal_dist ** 2

    positions = []
    velocities = []

    for gene in gene_key:
        z_all = data[data[0] == gene].values.tolist()[1:]

        position_arr = []
        velocity_arr = []

        # kalman loop
        for k in range(len(t)):
            z = z_all[k]
            x, P = predict(x, P, F, Q)
            x, P = update(x, P, z, R, H)

            dx = (t[k] - prev_x)
            F = np.array([[1.,dx],
                        [0.,1.]])
            prev_x = t[k]

            position_arr.append(x[0])
            velocity_arr.append(x[1])

        positions.append(position_arr)
        velocities.append(velocity_arr)
    
    velocities_matrix = np.array(velocities) # m genes, n time points
    positions_matrix = np.array(positions) # ^

    return positions_matrix, velocities_matrix, gene_key, t

def get_betas_for_timepoint(genes_position_vector, genes_velocity_vector, network_dict, gene_key, tf_key):
    coefficient_matrix = []
    beta_all_matrix = []

    for R in range(len(genes_velocity_vector)):
        gene = gene_key[R]

        # get beta_all 
        R_transcription =  genes_velocity_vector[R] + genes_position_vector[R] * decay_rate

        Kd_p = network_dict[gene_key[R]]["regulators"]["polymerase"]["Kd"]
        P_basal = Np / (Nns * Kd_p)
        R_basal = basal_rate

        beta_all = R_transcription / (P_basal * R_basal)
        beta_all_matrix.append(beta_all)

        # get coeffiecients
        tfs_info = [ tf for tf in network_dict[gene]["regulators"] ]

        coefficient_arr = [0 for i in len(tf_key)]

        for tf in tfs_info:
            N_tf = genes_position_vector[gene_key.index(str(tfs_info[tf].keys()))]
            P = N_tf / (N_tf + tf["kd"])
            coefficient_arr[tf_key.index(tf)] = P

        coefficient_matrix.append(coefficient_arr)

    # format matrices
    coefficient_matrix = np.array(coefficient_matrix)
    beta_all_matrix = np.array(beta_all_matrix)
    beta_all_matrix = np.reshape(beta_all_matrix, (-1, 1))

    # least squares fitting
    beta_vector = lsqr(coefficient_matrix, beta_all_matrix, iter_lim=1e6)
    beta_arr = beta_vector.values.flatten().tolist()

    return beta_arr

def main(prokode_dir, data_file, genome_loc, annotation_loc):
    # kalman filtering
    data = pd.read_csv(data_file)
    positions_matrix, velocities_matrix, gene_key, t = kalman_filtering(data)

    # create network.json
    network_loc = create_network_json_main(prokode_dir, genome_loc, annotation_loc)
    network_dict = json.load(open(network_loc, 'r'))

    # define arbitrary order of tfs
    tf_key = [ tf for tf in network_dict[i]["regulators"] for i in len(network_dict) ]
    tf_key = list(set(tf_key)) # filter repeats
    tf_key = tf_key.remove('polymerase') # filter out polymerase

    # get beta vals for each time point
    beta_collection = [] # m time points by n tfs
    for time_index in range(len(t)):
        beta_arr = get_betas_for_timepoint(positions_matrix[:, time_index], velocities_matrix[:, time_index], network_dict, gene_key, tf_key)
        beta_collection.append(beta_arr)
    beta_collection = np.array(beta_collection)

    compare(beta_collection)

main(prokode_dir, data_file, genome_loc, annotation_loc)