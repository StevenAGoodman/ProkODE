import json
import ast
import numpy as np
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

def create_network_json(prokode_dir, gene_arr, decay_loc, tfbs_loc, annotation_df):
    # file df dependencies: 
    decay_df = pd.read_csv(decay_loc) # names=['gene', 'decay_rate']
    tfbs_df = pd.read_csv(tfbs_loc).reset_index() # names=[,'tf','tg','kd', 'beta]

    # json start brackets
    output = {}

    # loop over every gene
    for tg in gene_arr:
        # query decay rates file
        tgdecay = decay_df.loc[decay_df['gene']==tg, 'decay_rate']

        # within each gene's brackets:
        print(tgdecay.to_list())
        syn = annotation_df.loc[annotation_df['geneid']==tg,'synonyms'].values[0]
        print(syn)
        arr = {"tg_decay":tgdecay.to_list()[0], "synonyms":ast.literal_eval(syn)}
        reg_arr = {}
        # within each gene's "regulators": array
        for _, tf_row in tfbs_df[tfbs_df['tg']==tg].iterrows():
            tf_row['tf'] = "polymerase" if tf_row['tf'] == "rpoD" else tf_row['tf']
            # beta = tf_row['beta']
            kdtf = tf_row['Kd']
            reg_arr[tf_row['tf']] = {"beta":"nan","kd_tf":kdtf}
        arr["regulators"] = reg_arr
        
        output[tg] = arr

    # write to network.json
    network_loc = prokode_dir + '/src/network.json'
    json.dump(output, open(network_loc, 'w'), indent=5)

    return network_loc

# def config_basal(Nns):
#     Np = 6000
#     Kdp =
#     return Np / (Nns * Kdp)

def network_main(prokode_dir, annotation_loc, tfbs_loc, decay_rates_loc):
    # configer basal / polymerase condition files
    # basal_loc =

    # create gene_arr needed for create_network_json()
    annotation_df = pd.read_csv(annotation_loc, delimiter='\t')
    gene_arr = list(set(annotation_df['geneid'].to_list()))

    # create network.json
    network_loc = create_network_json(prokode_dir, gene_arr, decay_rates_loc, tfbs_loc, annotation_df)

    return network_loc

# net = network_main('/workspaces/PROKODE-DOCKER', '/workspaces/PROKODE-DOCKER/src/inputs/annotation.csv', '/workspaces/PROKODE-DOCKER/src/tfbs.csv', '/workspaces/PROKODE-DOCKER/src/decay_rates.csv')