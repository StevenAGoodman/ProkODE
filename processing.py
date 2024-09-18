import pandas as pd
import numpy as np
from scipy.sparse.linalg import lsqr
import json
import statistics
from run import create_network_json_main
import cvxpy 

prokode_dir = 'C:/Users/cryst/LOFScreening/archive/PROKODE'
data_file = 'GSE90743_E14R025_raw_counts.txt'
genome_loc = 'C:/Users/cryst/LOFScreening/archive/PROKODE/src/inputs/genome.fasta'
annotation_loc = 'C:/Users/cryst/LOFScreening/archive/PROKODE/src/inputs/annotation.csv'
pfm_database_loc = 'C:/Users/cryst/LOFScreening/archive/PROKODE/src/preprocessing/pfmdb.txt'
CiiiDER_jar_loc = './CiiiDER_TFMs/CiiiDER.jar'
CiiiDER_thresh = 0.15

sensor_normal_dist = 10
basal_rate = 3
decay_rate = np.log(2)/300
Np = 6000
Kd_p = 0.1
Nns = 4600000

def find_key_nonrecursive(adict:dict, key):
    for k, v in adict.items():
        if key in k:
            return v
        d = v["synonyms"]
        if key in d:
            return v
    
def compare(beta_collection, tf_key):
    compare_json = []
    for tf_beta_index in range(len(tf_key)):
        tf = tf_key[tf_beta_index]
        tf_beta_vector = beta_collection[:, tf_beta_index]

        mean = statistics.fmean(tf_beta_vector)
        median_r = statistics.median(tf_beta_vector)
        range_r = [float(min(tf_beta_vector)), float(max(tf_beta_vector))]
        std_dev = statistics.stdev(tf_beta_vector)

        compare_json.append({"tf":tf,"mean":mean,"median":median_r,"range":range_r,"standard deviation":std_dev})
    
    json.dump(compare_json, open('results.json', 'w'), indent=3)

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
    gene_key = data.iloc[:, 0].values.flatten().tolist()

    # define t
    t = [5,10,25,45,75,120,210,330,1500,1560,1680]

    # init kalman inputs
    x = x_0 = np.array([[0],
                  [0]])
    prev_x = 0
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
        z_all = data.loc[data['Gene'] == gene].values[0][1:]
        position_arr = []
        velocity_arr = []

        # kalman loop
        for k in range(len(t)):
            z = z_all[k]
            x, P = predict(x, P, F, Q)
            x, P = update(P, H, R, z, x)

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

    print('kalman complete...')
    return positions_matrix, velocities_matrix, gene_key, t

def get_betas_for_timepoint(genes_position_vector, genes_velocity_vector, network_dict, gene_key, tf_key):
    coefficient_matrix = []
    beta_all_matrix = []

    for R in range(len(genes_velocity_vector)):
        gene = gene_key[R]

        # get beta_all 
        R_transcription =  genes_velocity_vector[R] + genes_position_vector[R] * decay_rate

        # Kd_p = network_dict[gene_key[R]]["regulators"]["polymerase"]["Kd"]
        P_basal = Np / (Nns * Kd_p)
        R_basal = basal_rate

        beta_all = R_transcription / (P_basal * R_basal)
        beta_all_matrix.append(beta_all)

        # get coeffiecients
        k = find_key_nonrecursive(network_dict, gene)
        try:
            tfs_info = [ tf for tf in k["regulators"].keys() ]
        except:
            print(gene, 'FAILURE!!!!!')
            tfs_info = []
        coefficient_arr = [0 for i in range(len(tf_key))]

        for tf in tfs_info:
            N_tf = genes_position_vector[gene_key.index(tf)]
            P = N_tf / (N_tf + k["regulators"][tf]["kd_tf"])
            coefficient_arr[tf_key.index(tf)] = float(P)
        
        coefficient_matrix.append(coefficient_arr)

    # format matrices
    coefficient_matrix = np.array(coefficient_matrix)
    beta_all_matrix = np.array(beta_all_matrix)
    beta_all_matrix = np.reshape(beta_all_matrix, (-1, 1))

    # least squares fitting
    print('\t...solving least squares problem')
    X = cvxpy.Variable((coefficient_matrix.shape[1],1))
    constraints = [X >= 0]

    product = coefficient_matrix @ X
    cost = cvxpy.sum(product, axis=1, keepdims=True) - beta_all_matrix
    problem = cvxpy.Problem(cvxpy.Minimize(cvxpy.norm(cost)), constraints)
    problem.solve()

    # beta_vector = solve_ls(coefficient_matrix, beta_all_matrix, lb=lower_bound, solver='osqp', sparse_conversion=True, time_limit=60)
    # print(beta_vector)
    beta_arr = X.value.transpose().tolist()[0]
    print(beta_arr)

    return beta_arr

def main(prokode_dir, data_file, genome_loc, annotation_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, prebuilt=False):
    print("=====================================\nProkODE: GetBetas has begun. Running...\n=====================================\n\n")

    # create network.json
    print('creating network.json...')
    if prebuilt:
        network_loc = prokode_dir + '/src/network.json'
    else:
        network_loc = create_network_json_main(prokode_dir, genome_loc, annotation_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, False)# 
    network_dict: dict = json.load(open(network_loc, 'r'))
    

    # define arbitrary order of tfs
    tf_key = [ tf for val in network_dict.values() for tf in val["regulators"].keys() ]
    tf_key = list(set(tf_key)) # filter repeats

    # kalman filtering
    print('kalman filtering...')
    data = pd.read_csv(data_file, delimiter='\t')
    positions_matrix, velocities_matrix, gene_key, t = kalman_filtering(data)

        
    # get beta vals for each time point
    print('getting beta collection...')
    beta_collection = [] # m time points by n tfs
    for time_index in range(len(t)):
        print('\t...getting betas for timepoint')
        beta_arr = get_betas_for_timepoint(positions_matrix[:, time_index], velocities_matrix[:, time_index], network_dict, gene_key, tf_key)
        beta_collection.append(beta_arr)
    print('\t...merging timepoints')
    beta_collection = np.array(beta_collection)

    print('analyzing and writing results file...')

    compare(beta_collection, tf_key)

main(prokode_dir, data_file, genome_loc, annotation_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, False)