### UNDER DEVELOPMENT ###
import math
import json
import linecache
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from maths import *
from decay.processing import decay_main

def search_json(json_loc, key, gene_key):
    index = gene_key.index(key)
    line = linecache.getline(json_loc, index)
    line = line[line.find('{'):-2]
    gene_dict = json.loads(line)
    return gene_dict

def get_decay_dicts():

    return {}, {}

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
