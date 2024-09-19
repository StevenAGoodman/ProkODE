### UNDER DEVELOPMENT ###
import json
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from maths import * 

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

def update_amnts(gene_key:list, prev_mRNA_amnts:list, prev_protein_amnts:list, dt:float):
    mRNA_amnts = []
    protein_amnts = []

    for n in range(len(gene_key)):
        gene = gene_key[n]

        # mRNA calculation
        mRNA_creation_rate = transcription_rate(gene, prev_protein_amnts)
        mRNA_decay_rate = decay_rate("mRNA", gene, prev_protein_amnts)
        total_mRNA_rate = mRNA_creation_rate - mRNA_decay_rate * prev_mRNA_amnts[n]
        mRNA_amnt = prev_mRNA_amnts[n] + total_mRNA_rate * dt

        # protein calculation
        prot_creation_rate = translation_rate(gene, prev_protein_amnts)
        prot_decay_rate = decay_rate("protein", gene, prev_protein_amnts)
        total_prot_rate = prot_creation_rate - prot_decay_rate * prev_protein_amnts[n]
        protein_amnt = prev_protein_amnts[n] + total_prot_rate * dt

        mRNA_amnts.append(mRNA_amnt)
        protein_amnts.append(protein_amnt)

    return mRNA_amnts, protein_amnts

def get_plotting_data(gene_key, init_mRNA_amnts, init_protein_amnts, max_time, dt):
    x_coords = [0]
    y_coords = [init_protein_amnts]
    
    time = 0
    mRNA_amnts = init_mRNA_amnts
    protein_amnts = init_protein_amnts

    for i in math.floor(max_time / dt):
        time += dt
        x_coords.append(time)
        
        mRNA_amnts, protein_amnts = update_amnts(gene_key, mRNA_amnts, protein_amnts, dt)

        y_coords.append(protein_amnts)

    return x_coords, y_coords

def plot_system(gene_key, init_mRNA_amnts, init_protein_amnts, dt, max_time):
    x_coords, y_coords_tuple = get_plotting_data(gene_key, init_mRNA_amnts, init_protein_amnts, max_time, dt)

    for i in range(10):
        plt.plot(t,s[:,i],random.choice(['r-','g-','b-']), linewidth=2.0)
    # plt.plot(t,s[:,2], 'g-',linewidth=2.0)
    plt.xlabel("t")
    plt.ylabel("S[N,C]")
    plt.legend(["N","C",'d'])
    plt.show()
