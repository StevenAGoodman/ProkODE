import json
import ast
import re
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

def create_network_json(prokode_dir, tfbs_loc, annotation_df, operons_df):
    # json start brackets
    output = {}
    gene_key = []

    # loop over every gene
    for i, _ in operons_df.iterrows():
        operon_gene_arr:str = operons_df.loc[i, 'geneids']
        operon_gene_arr:list = json.loads(operon_gene_arr)

        # length of multi-gene mrna from operon
        transcript_len = int(operons_df.loc[i, 'end']) - int(operons_df.loc[i, 'start'])

        for gene in operon_gene_arr:
            try: # sometimes geneids from the operon.tsv are not compatible with geneids from annotation.tsv
                # length of individual gene mrna
                mRNA_len = int(annotation_df.loc[annotation_df['geneid']==gene, 'end'].tolist()[0]) - int(annotation_df.loc[annotation_df['geneid']==gene, 'start'].tolist()[0])

                # within each gene's brackets:
                syn = annotation_df.loc[annotation_df['geneid']==gene,'synonyms'].tolist()[0]
                arr = {"synonyms":ast.literal_eval(syn),"transcript length":transcript_len,"mRNA length":mRNA_len}
                reg_arr = {}
                arr["regulators"] = reg_arr
                
                output[gene] = arr
                gene_key.append(gene)
            except:
                print('failed: ', gene)
                None
    
    with open(tfbs_loc, 'r') as tfbs_file:
        for _, line in enumerate(tfbs_file):
            # skip blank lines
            if re.search(r'.+\d', line) == None:
                continue

            line = line.split(',')
            operon = int(line[1])
            operon_genes:str = operons_df.loc[operons_df['operonid']==operon, 'geneids'].tolist()[0]
            operon_genes:list = json.loads(operon_genes)

            tf = line[0]
            kdtf = float(line[2].replace('\n',''))
            # beta = line[3]
            
            try:
                for tg in operon_genes:
                    output[tg]["regulators"][tf] = {"beta":"NaN","score":kdtf}
            except:
                continue

    # write to network.json
    print('\t\twriting to network.json')
    network_loc = prokode_dir + '/src/network.json'
    with open(network_loc, 'a') as network_file:
        network_file.write('{\n')
        for gene, gene_elmnt in output.items():
            if gene == list(output.keys())[-1]: # make sure theres no comma after the last item
                network_file.write(f'"{gene}"' + ':' + json.dumps(gene_elmnt) + '\n')
            else:
                network_file.write(f'"{gene}"' + ':' + json.dumps(gene_elmnt) +',\n')
        network_file.write('}')

    return network_loc, gene_key

# def create_network_json(prokode_dir, gene_arr, decay_loc, tfbs_loc, annotation_df):
#     # read file dependencies: 
#     decay_df = pd.read_csv(decay_loc) # names=['gene', 'decay_rate']
    
#     # read tfbs.csv and append to network.json
#     print('\t\twriting to network.json')
#     network_loc = prokode_dir + '/src/network.json'
#     with open(tfbs_loc, 'r') as tfbs_file: # names=[,'tf','tg','kd', 'beta]
#         with open(network_loc, 'a') as network_file:
#             # json start bracket
#             network_file.write('{\n')

#             # loop over every gene
#             for tg in gene_arr:
#                 # query decay rates file
#                 tgdecay = decay_df.loc[decay_df['gene']==tg, 'decay_rate']

#                 # within each gene's brackets:
#                 syn = annotation_df.loc[annotation_df['geneid']==tg,'synonyms'].tolist[0]
#                 arr = {"tg":tg, "synonyms":ast.literal_eval(syn), "tg_decay":tgdecay.to_list()[0]}
#                 reg_arr = {}
#                 # within each gene's "regulators": array
#                 for _, tf_row in tfbs_df[tfbs_df['tg']==tg].iterrows():
#                     tf_row['tf'] = "polymerase" if tf_row['tf'] == "rpoD" else tf_row['tf']
#                     # beta = tf_row['beta']
#                     kdtf = tf_row['Kd']
#                     reg_arr[tf_row['tf']] = {"beta":"nan","kd_tf":kdtf}
#                 arr["regulators"] = reg_arr
                
#                 network_file.write(f'{json.loads(arr)},\n')
            
#             # json end brackets
#             network_file.write('}')

#     return network_loc

# def config_basal(Nns):
#     Np = 6000
#     Kdp =
#     return Np / (Nns * Kdp)

def network_main(prokode_dir, annotation_loc, operons_loc, tfbs_loc):
    # configer basal / polymerase condition files
    # basal_loc =

    # create gene_arr needed for create_network_json()
    annotation_df = pd.read_csv(annotation_loc, delimiter='\t')
    operons_df = pd.read_csv(operons_loc, sep='\t')

    # create network.json
    network_loc, gene_key = create_network_json(prokode_dir, tfbs_loc, annotation_df, operons_df)

    return network_loc, gene_key

# net = network_main('/workspaces/PROKODE-DOCKER', '/workspaces/PROKODE-DOCKER/src/inputs/annotation.csv', '/workspaces/PROKODE-DOCKER/src/tfbs.csv', '/workspaces/PROKODE-DOCKER/src/decay_rates.csv')