import cvxpy.atoms
import cvxpy.atoms.elementwise
import cvxpy.atoms.elementwise.log
import pandas as pd
import numpy as np
import json
import statistics
from matplotlib import pyplot as plt
# from src.run import create_network_json_main
from src.plot_system import *
import cvxpy 

reset = True
prokode_dir = '/workspaces/PROKODE'
data_file = 'GSE90743_E14R025_raw_counts.txt'
genome_loc = prokode_dir + '/src/inputs/genome.fasta'
annotation_loc = prokode_dir + '/src/inputs/annotation.tsv'
operons_loc = prokode_dir + '/src/inputs/operons.tsv'
pfm_database_loc = prokode_dir + '/src/preprocessing/pfmdb.txt'
CiiiDER_jar_loc = './CiiiDER_TFMs/CiiiDER.jar'
CiiiDER_thresh = 0.3

# global jazz
sensor_normal_dist = 10
basal_rate = 3
decay_rate = np.log(2)/300
Np = 6000
Kd_rnap = 0.1
Nns = 4600000
genome_len = 4.64e6

# constants
temperature = 200 # kelvin
elongation_rate = 60 # nt/s
peptide_rate = 20 # aa/s
Kd_ribo_mrna = 0.1 # WHAT IS THE BINDING AFFINITY OF RIBO TO DALGARNO SEQ????
len_taken_by_rnap = 30 # nt
len_taken_by_ribo = 30 # aa

def update_mean(count, existing_mean, new_value):
    count += 1
    delta = new_value - existing_mean
    mean = existing_mean + (delta / count)
    return (count, mean)

def find_key_nonrecursive(adict:dict, key):
    for k, v in adict.items():
        if key in k:
            return v
        d = v["synonyms"]
        if key in d:
            return v
    
def compare(t, beta_collection, tf_key):
    compare_json = []
    count = 0
    running_mean_std_dev = 0
    for tf_beta_index in range(len(tf_key)):
        tf = tf_key[tf_beta_index]
        tf_beta_vector = beta_collection[:, tf_beta_index]

        mean = statistics.fmean(tf_beta_vector)
        median_r = statistics.median(tf_beta_vector)
        range_r = [float(min(tf_beta_vector)), float(max(tf_beta_vector))]
        std_dev = statistics.stdev(tf_beta_vector)

        count, running_mean_std_dev = update_mean(count, running_mean_std_dev, std_dev)

        compare_json.append({"tf":tf,"mean":mean,"median":median_r,"range":range_r,"log standard deviation":np.log(std_dev)})
    
    compare_json.append({"tf":"all","mean standard deviation":running_mean_std_dev,"mean log standard deviation":np.log(running_mean_std_dev)})
    
    json.dump(compare_json, open('results.json', 'w'), indent=3)

    y = []
    for i in range(len(beta_collection)):
        std_dev = statistics.stdev(beta_collection[i])
        y.append(std_dev)
    
    plt.plot(t, y, '-b')
    plt.show()

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
    t = [60 * min for min in t] # to seconds

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

def get_betas_for_timepoint(prev_protein_amnts, genes_position_vector, genes_velocity_vector, N_rnap, N_ribo, network_loc, network_gene_key, mRNA_decay_dict, protein_decay_dict, gene_key:list, tf_key, dt):
    coefficient_matrix = []
    beta_all_matrix = []

    protein_amnts = []

    for R in range(len(genes_velocity_vector)):
        gene = gene_key[R]

        # get coeffiecients
        try:
            k = search_json(network_loc, gene, network_gene_key)
            tfs_info = [ tf for tf in k["regulators"].keys() ]
        except:
            print(gene, 'FAILURE!!!!!')
            tfs_info = []
        coefficient_arr = [0 for i in range(len(tf_key))]

        for tf in tfs_info:
            N_tf = prev_protein_amnts[gene_key.index(tf)]
            P = N_tf / (N_tf + score_to_K(k["regulators"][tf]["score"]))

            if float(P) == 0: # opportunity to replace with psuedocounts
                P += 0.000001

            # print(k["regulators"][tf]["score"], score_to_K(k["regulators"][tf]["score"]), P)
            coefficient_arr[tf_key.index(tf)] = float(P) # ~multiplicative method: np.log10(float(P))
        # get mRNA decay
        mRNA_decay_rate = RNA_decay_rate(gene, genes_position_vector, genes_position_vector[R], mRNA_decay_dict, gene_key)

        # get beta_all 
        R_transcription = genes_velocity_vector[R] + genes_position_vector[R] * mRNA_decay_rate

        # calculate max transcription rate
        # transcript_len = gene_info_dict["transcript length"]
        max_txn_rate = 1 / (len_taken_by_rnap / elongation_rate) # transcripts/s

        # calculate basal rnap binding prob
        P_rnap_basal = N_rnap / (genome_len * (N_rnap + Kd_rnap))

        beta_all = R_transcription / (P_rnap_basal * max_txn_rate)
        beta_all_matrix.append(np.log(beta_all)) # ~multiplicative method: np.log10(beta_all) - sum(coefficient_arr)

        # get protein rates
        if k != []:
            prot_creation_rate = translation_rate(gene, prev_protein_amnts, N_ribo, k)
            prot_decay_rate = protein_decay_rate(gene, prev_protein_amnts, protein_decay_dict, gene_key)
            total_prot_rate = prot_creation_rate * genes_position_vector[R] - prot_decay_rate * prev_protein_amnts[R]
            protein_amnt = prev_protein_amnts[R] + total_prot_rate * dt
        else:
            protein_amnt = 0

        protein_amnts.append(protein_amnt)

        coefficient_matrix.append(coefficient_arr) # mult: [ 1 for i in range(len(tf_key))]

    # format matrices
    coefficient_matrix = np.array(coefficient_matrix)
    beta_all_matrix = np.array(beta_all_matrix)
    beta_all_matrix = np.reshape(beta_all_matrix, (-1, 1))
    print(coefficient_matrix, '\n', beta_all_matrix)

    # least squares fitting
    print('\t...solving least squares problem')
    # X = cvxpy.Variable((coefficient_matrix.shape[1],1))
    # constraints = [X >= 0]

    # cost = cvxpy.log(product) - beta_all_matrix
    # problem = cvxpy.Problem(cvxpy.Minimize(cvxpy.norm(cost, p=2)), constraints)
    # problem.solve(verbose=True)
    # Define variables
    coefficient_matrix = np.random.rand(100, 10)
    beta_all_matrix = np.random.rand(100, 1)

    log_beta = cvxpy.Variable((coefficient_matrix.shape[1], 1))  # log of beta values
    constraints = []

    # Convert to a linear form: log(P) + log(beta) = log(beta_all)
    log_product = cvxpy.log(coefficient_matrix) @ log_beta
    cost = cvxpy.norm(log_product - cvxpy.log(beta_all_matrix), 2) + 1e-6 * cvxpy.norm(log_beta)

    # Formulate the problem
    problem = cvxpy.Problem(cvxpy.Minimize(cost), constraints)
    problem.solve( verbose=True)

    # Obtain original beta values
    beta_values = np.exp(log_beta.value)
    print(beta_values)

    # beta_vector = solve_ls(coefficient_matrix, beta_all_matrix, lb=lower_bound, solver='osqp', sparse_conversion=True, time_limit=60)
    # print(beta_vector)
    beta_arr = X.value.transpose().tolist()[0]
    print(beta_arr)

    return beta_arr, protein_amnts

def main(prokode_dir, data_file, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, prebuilt=False, reset=True):
    print("=====================================\nProkODE: GetBetas has begun. Running...\n=====================================\n")

    # create network.json if not already built
    print('creating network.json...')
    if prebuilt:
        network_loc = prokode_dir + '/src/network.json'
    else:
        None
        # network_loc = create_network_json_main(prokode_dir, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, False, reset)# 
        # define arbitrary order of tfs
    tf_key = [ tf for val in json.load(open(network_loc, 'r')).values() for tf in val["regulators"].keys() ]
    tf_key = list(set(tf_key)) # filter repeats

    network_gene_key = []
    with open(network_loc, 'r') as network_f:
        for _, line in enumerate(network_f):
            if line == '{\n' or line == '}':
                continue
            else:
                gene = line[line.find('"')+1:line.find('":{')]
                network_gene_key.append(gene)

    # kalman filtering
    print('kalman filtering...')
    data = pd.read_csv(data_file, delimiter='\t')
    positions_matrix, velocities_matrix, gene_key, t = kalman_filtering(data)

        
    # get beta vals for each time point
    print('getting beta collection...')
    beta_collection = [] # m time points by n tfs
    protein_amnts = positions_matrix[:, 0]
    mRNA_decay_dict, protein_decay_dict = get_decay_dicts()
    N_rnap = 10000
    N_ribo = 15000
    time = 0

    for time_index in range(len(t)):
        dt = t[time_index] - time
        time = t[time_index]

        print('\t...getting betas for timepoint')
        beta_arr, protein_amnts = get_betas_for_timepoint(protein_amnts, positions_matrix[:, time_index], velocities_matrix[:, time_index], N_rnap, N_ribo, network_loc, network_gene_key, mRNA_decay_dict, protein_decay_dict, gene_key, tf_key, dt)
        beta_collection.append(beta_arr)

    print('\t...merging timepoints')
    beta_collection = np.array(beta_collection)

    print('analyzing and writing results file...')

    compare(t, beta_collection, tf_key)

main(prokode_dir, data_file, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, True, reset)