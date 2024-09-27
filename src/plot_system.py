### UNDER DEVELOPMENT ###
import json
import linecache
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from src.maths import * 
from src.decay.processing import decay_main

# constants
temperature = 298 # kelvin
elongation_rate = 60 # nt/s
peptide_rate = 20 # aa/s
Kd_ribo_mrna = 0.1 # WHAT IS THE BINDING AFFINITY OF RIBO TO DALGARNO SEQ????
len_taken_by_rnap = 30 # nt
len_taken_by_ribo = 30 # aa

def search_json(json_loc, key, gene_key):
    index = gene_key.index(key)
    line = linecache.getline(json_loc, index)
    line = line[line.find('{'):-2]
    gene_dict = json.loads(line)
    return gene_dict

def score_to_K(score):
    delta_G = -1 * score * 1000
    return 1 / np.exp(delta_G / (1.98722 * temperature))

def get_decay_dicts():

    return {}, {}

def transcription_rate(gene, protein_amnts, gene_key, N_rnap, Kd_rnap, gene_info_dict, genome_len):
    genome_len = 4.5e6

    # calculate beta_all
    beta_all = 0.0
    for tf, tf_info in gene_info_dict["regulators"].items():
        N_tf = protein_amnts[gene_key.index(tf)]
        P = N_tf / (N_tf + score_to_K(tf_info["score"]))
        beta = tf_info["beta"]
        beta_all += float(beta) * float(P)

    # calculate max transcription rate
    transcript_len = gene_info_dict["transcript length"]
    max_txn_rate = 1 / (len_taken_by_rnap / elongation_rate) # transcripts/s

    # calculate basal rnap binding prob
    P_rnap_basal = N_rnap / (genome_len * (N_rnap + Kd_rnap))

    return beta_all * P_rnap_basal * max_txn_rate

def translation_rate(gene, protein_amnts, N_ribo, gene_info_dict):
    if gene_info_dict == None:
        return 0
    # calculate max translation rate
    mRNA_len = gene_info_dict["mRNA length"]
    max_translation_rate = 1 / (len_taken_by_ribo / peptide_rate)

    # calculate ribo binding prob
    P_ribo_bound = N_ribo / (N_ribo + Kd_ribo_mrna)

    return P_ribo_bound * max_translation_rate

def RNA_decay_rate(gene, prev_total_mRNA_amnt, protein_amnts, decay_dict, gene_key):
    # calculation of natural decay component
    natural_mRNA_decay_rate = 0

    # get degrad_prots of gene
    decay_info = decay_dict#[gene] # should be a dict of 

    rate_of_mRNA_cleavage = 0
    for degrad_prot, K in decay_info.items():
        N_dp = protein_amnts[gene_key.index(degrad_prot)]
        # K_1 is the reaction rate from [Enzyme] + [mRNA] -> [Enzyme-mRNA] (ie, K_1[Enzyme][mRNA] = [Enzyme-mRNA] create / time ) ... k_1 is in 1 / (mols * time)
        rate_of_mRNA_cleavage += K * prev_total_mRNA_amnt * N_dp / (1 + K * prev_total_mRNA_amnt)
        # what is the relationship of multiple degrading proteins? it should be additive, right?
        # how can you identify a protien as degrading?
        # what is the relationship of degrading prots component and decay (ie half life) component

    return rate_of_mRNA_cleavage + natural_mRNA_decay_rate

def protein_decay_rate(gene, protein_amnts, decay_dict, gene_key):
    # total_half_life = 0

    # # N end rule
    # total_half_life += 1 / decay_dict["N end"][gene]

    # # decay from misfolded proteins
    # total_half_life += 1 / decay_dict["misfold"][gene]()

    # return 1 / total_half_life
    return 0

def update_amnts(gene_key:list, prev_mRNA_amnts:list, prev_protein_amnts:list, N_rnap, N_ribo, mRNA_decay_dict, protein_decay_dict, network_loc, dt:float):
    mRNA_amnts = []
    protein_amnts = []

    for n in range(len(gene_key)):
        gene = gene_key[n]
        gene_info_dict = linecache.getline(network_loc, n)
        gene_info_dict = json.loads(gene_info_dict)

        # mRNA calculation
        mRNA_creation_rate = transcription_rate(gene, prev_protein_amnts, gene_key, N_rnap, 0.1, gene_info_dict, genome_len=0)
        mRNA_decay_rate = RNA_decay_rate(gene, prev_mRNA_amnts, prev_protein_amnts, mRNA_decay_dict, gene_key)    
        total_mRNA_rate = mRNA_creation_rate - mRNA_decay_rate
        mRNA_amnt = prev_mRNA_amnts[n] + total_mRNA_rate * dt

        # protein calculation
        prot_creation_rate = translation_rate(gene, prev_protein_amnts, N_ribo, gene_info_dict)
        prot_decay_rate = protein_decay_rate(gene, prev_protein_amnts, protein_decay_dict, gene_key)
        total_prot_rate = prot_creation_rate * prev_mRNA_amnts[n] - prot_decay_rate * prev_protein_amnts[n]
        protein_amnt = prev_protein_amnts[n] + total_prot_rate * dt

        mRNA_amnts.append(mRNA_amnt)
        protein_amnts.append(protein_amnt)

    # update rnap and ribosome populations

    return mRNA_amnts, protein_amnts, N_rnap, N_ribo

def get_plotting_data(init_mRNA_amnts, init_protein_amnts, network_loc, gene_key, mRNA_decay_dict, protein_decay_dict, max_time, dt):
    x_coords = [0]
    y_coords = [init_protein_amnts]

    # initialize variables
    time = 0
    mRNA_amnts = init_mRNA_amnts
    protein_amnts = init_protein_amnts
    N_rnap = 10000
    N_ribo = 15000

    for i in range(math.floor(max_time / dt)):
        time += dt
        x_coords.append(time)
        
        mRNA_amnts, protein_amnts, N_rnap, N_ribo = update_amnts(gene_key, mRNA_amnts, protein_amnts, N_rnap, N_ribo, mRNA_decay_dict, protein_decay_dict, network_loc, dt)

        y_coords.append(protein_amnts)

    return x_coords, y_coords

def plot_system(prokode_dir, network_loc, mRNA_decay_loc, protein_decay_loc, gene_key, dt, max_time):
    print('== plotting begun ==')
    # mRNA_decay_loc, protein_decay_loc = decay_main(prokode_dir, gene_key)


    # default to no initial mRNA
    # try: read file except:
    init_mRNA_amnts = [0] * len(gene_key)

    # read input init protein file
    init_protein_amnts = [np.random.randint(0,1000) for i in range(len(gene_key))]

    # config decay_dicts\
    mRNA_decay_dict, protein_decay_dict = get_decay_dicts()

    x_coords, y_coords_tuple = get_plotting_data(init_mRNA_amnts, init_protein_amnts, network_loc, gene_key, mRNA_decay_dict, protein_decay_dict, max_time, dt)

    for i in range(10):
        plt.plot(x_coords,[gene_exps[i] for gene_exps in y_coords_tuple], random.choice(['r-','g-','b-']), linewidth=2.0)
    # plt.plot(t,s[:,2], 'g-',linewidth=2.0)
    plt.xlabel("t")
    plt.ylabel("S[N,C]")
    plt.legend(["N","C",'d'])
    plt.show()

# prokode_dir = 'C:/Users/cryst/LOFScreening/archive/PROKODE'
# plot_system(prokode_dir, prokode_dir + '/src/network.json',1,1,1,100)

# network.json formatting:
# {
#     "tg1":[
#         "tgdecay":float,
#         "regulators":[
#             "tf1": [
#                 "beta":float,
#                 "kd_tf":float,
#             ]
#             "tf2":[

#             ]
#             ...
#         ]
#     ],
#     "tg2":[

#     ],
#     ...
# } 

# decay_dict formating:
# {
#     "deg_enzyme1":{
#         "targe protein":lambda deg_enzyme_amnt, target_prot_amnt : xyz
#         "target protein":lambda:
#     }
# }