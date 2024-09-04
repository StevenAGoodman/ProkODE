import json
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd

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

def config_network_json(decay_file, tfbs_file):
    # needed inputs: lists of tgs, gene - decay rates ref file,

    # implicit file dependencies: 
    decay_df = pd.read_csv(decay_file) # names=['gene', 'prot decay rate']
    topology_df = pd.read_csv(tfbs_file) # names=['tf','tg','kd', 'beta]
    topology_df = topology_df.reset_index()

    tg_arr = list(set(topology_df['tg'].to_list()))
    output = {}

    for tg in tg_arr:
        print(tg)
        tgdecay = decay_df.at[decay_df['gene']==tg, 'decay']# query decay rates file
        
        arr = {"tg_decay":tgdecay,}
        reg_arr = {}
        for _, tf_row in topology_df[topology_df['tg']==tg].iterrows():
            beta = tf_row['beta']
            kdtf = tf_row['kd']
            reg_arr[tf_row['tf']] = {"beta":beta,"kd_tf":kdtf}
        arr["regulators"] = reg_arr
        print(arr)
        output[tg] = arr

    json.dump(output, open('network.json', 'w'))

    return tg_arr

# def get_basal(Nns):
    # Np = 
    # Kdp =
    return 0. # Np / (Nns * Kdp)

def plot_system(params):
    # check if params len mathces nodes

    def f(s,t):
        diffeqs = []
        prev_ids = [] # corresponds to prev_amnts; use as reference for index of id

        for tg in (params):
            param = params[i+1]
            R_trans = (prev_amnts[i] / (prev_amnts[i] + param[0])) * (param[1] * E_basal) + (1 - (prev_amnts[i] / (prev_amnts[i] + param[0]))) * E_basal
            d = R_trans - param[2] * s[i+1]
            prev_amnts.append(d)
                
        return diffeqs

    t = np.linspace(0,10000)
    s0=[0,0]

    s = odeint(f,s0,t)

    plt.plot(t,s[:,0],'r-', linewidth=2.0)
    plt.plot(t,s[:,1],'b-', linewidth=2.0)
    # plt.plot(t,s[:,2], 'g-',linewidth=2.0)
    plt.xlabel("t")
    plt.ylabel("S[N,C]")
    plt.legend(["N","C",'d'])
    plt.show()

tg_arr = config_network_json('decay_rates.csv', 'tfbs.csv')
params = json.load(open('network.json', 'r'))
s0 = [np.random.randint(0,100) for i in range(len(tg_arr))]
E_basal = 0.00434782608696
# plot_system(params, s0)