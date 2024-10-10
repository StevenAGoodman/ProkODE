import json
import ast
import re
import math
import linecache
import numpy as np
import pandas as pd
from Bio.Seq import Seq
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

n_end_rule = {
    "R": 2 * 60,
    "K": 2 * 60,
    "F": 2 * 60,
    "L": 2 * 60,
    "W": 2 * 60,
    "Y": 2 * 60,
    "H": 10 * 3600,
    "I": 10 * 3600,
    "D": 10 * 3600,
    "E": 10 * 3600,
    "N": 10 * 3600,
    "Q": 10 * 3600,
    "C": 10 * 3600,
    "A": 10 * 3600,
    "S": 10 * 3600,
    "T": 10 * 3600,
    "G": 10 * 3600,
    "V": 10 * 3600,
    "P": 10 * 3600,
    "M": 10 * 3600    
}

def get_N_end(gene, annotation_df, genome_loc):
    end_loc = int(annotation_df.loc[annotation_df['geneid']==gene, "end"])
    print(end_loc)
    line = linecache.getline(genome_loc, 1 + math.floor(end_loc / 81)).replace('\n','')
    print(line)
    rem = end_loc % 81
    print(rem)
    if rem > 3:
        n_residue = Seq(line[rem - 3: rem])
        print(n_residue)
        n_residue = str(n_residue.translate())
    else:
        n_residue = linecache.getline(genome_loc, math.floor(end_loc / 81)).replace('\n','')[-3 + rem:] + line[:rem]
        print(n_residue)
        n_residue = str(Seq(n_residue).translate())
    
    if n_residue == "*":
        end_loc
        if rem > 3:
            n_residue = Seq(line[rem - 3: rem])
            print(n_residue)
            n_residue = str(n_residue.translate())
        else:
            n_residue = linecache.getline(genome_loc, math.floor(end_loc / 81)).replace('\n','')[-3 + rem:] + line[:rem]
            print(n_residue)
            n_residue = str(Seq(n_residue).translate())

        

    return n_end_rule[n_residue]

def create_network_json(prokode_dir, tfbs_loc, annotation_df, operons_df, genome_loc, floating_genes):
    # json start brackets
    output = {}
    gene_key = []

    mrna_decay_prots = {"rne":0, "rnc":0, "rbn":0, "rng":0} # binding afinities to mrna
    prot_decay_prots = {"clpX":0, "clpP":0, "lon":0} # binding affinities to proteins
    # prot_decay_prots = {"clpX":{"kd":-, "targets":["",""]},}

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

                # gene decays
                mRNA_decay_arr = {}
                for i in range(len(mrna_decay_prots)):
                    pair = list(mrna_decay_prots.items())[i]
                    if pair[0] in list(annotation_df["geneid"]):
                        mRNA_decay_arr[pair[0]] = pair[1]
                protein_decay_arr = {}
                for i in range(len(prot_decay_prots)):
                    pair = list(prot_decay_prots.items())[i]
                    if pair[0] in list(annotation_df["geneid"]):
                        protein_decay_arr[pair[0]] = pair[1]

                n_end_rule = get_N_end(gene, annotation_df, genome_loc)

                # within each gene's brackets:
                syn = annotation_df.loc[annotation_df['geneid']==gene,'synonyms'].tolist()[0]
                arr = {"synonyms":ast.literal_eval(syn),"transcript length":transcript_len,"mRNA length":mRNA_len, "mRNA decay": mRNA_decay_arr, "protein decay": protein_decay_arr, "N end rule": n_end_rule}
                reg_arr = {"polymerase":8.1}
                arr["regulators"] = reg_arr

                output[gene] = arr
                gene_key.append(gene)
            except:
                print('\t\tfailed:', gene)
                None
    for gene in floating_genes:
        try: # sometimes geneids from the operon.tsv are not compatible with geneids from annotation.tsv
            # length of individual gene mrna
            mRNA_len = int(annotation_df.loc[annotation_df['geneid']==gene, 'end'].tolist()[0]) - int(annotation_df.loc[annotation_df['geneid']==gene, 'start'].tolist()[0])
            transcript_len = mRNA_len

            # gene decays
            mRNA_decay_arr = {}
            for i in len(mrna_decay_prots):
                pair = list(mrna_decay_prots.items())[i]
                if pair[0] in list(annotation_df["geneid"]):
                    mRNA_decay_arr[pair[0]] = pair[1]
            protein_decay_arr = {}
            for i in range(len(prot_decay_prots)):
                pair = list(prot_decay_prots.items())[i]
                if pair[0] in list(annotation_df["geneid"]):
                    protein_decay_arr[pair[0]] = pair[1]

            # within each gene's brackets:
            syn = annotation_df.loc[annotation_df['geneid']==gene,'synonyms'].tolist()[0]
            arr = {"synonyms":ast.literal_eval(syn),"transcript length":transcript_len,"mRNA length":mRNA_len, "mRNA decay": mRNA_decay_arr, "protein decay": protein_decay_arr}
            reg_arr = {"polymerase":8.1}
            arr["regulators"] = reg_arr

            output[gene] = arr
            gene_key.append(gene)
        except:
            print('\t\tfailed:', gene)
            None

    with open(tfbs_loc, 'r') as tfbs_file:
        for _, line in enumerate(tfbs_file):
            # skip blank lines
            if re.search(r'.+\d', line) == None:
                continue

            line = line.split(',')
            if len(line) < 2:
                print("line err:", line)
                continue

            try:
                operon = int(line[1])
                operon_genes:str = operons_df.loc[operons_df['operonid']==operon, 'geneids'].tolist()[0]
                operon_genes:list = json.loads(operon_genes)

                tf = line[0]
                kdtf = float(line[2].replace('\n',''))
                # beta = line[3]

                try:
                    for tg in operon_genes:
                        output[tg]["regulators"][tf] = {"beta":"NaN","delta G":kdtf}
                except:
                    continue
            except:
                tg = line[1]
                try:
                    output[tg]["regulators"][tf] = {"beta":"NaN","delta G":kdtf}
                except:
                    continue

    # write to network.json
    print('\twriting to network.json')
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


def network_main(prokode_dir, genome_loc, annotation_loc, operons_loc, tfbs_loc, floating_genes):
    # configer basal / polymerase condition files
    # basal_loc =

    # create gene_arr needed for create_network_json()
    annotation_df = pd.read_csv(annotation_loc, delimiter='\t')
    operons_df = pd.read_csv(operons_loc, sep='\t')

    # create network.json
    network_loc, gene_key = create_network_json(prokode_dir, tfbs_loc, annotation_df, operons_df, genome_loc, floating_genes)

    return network_loc, gene_key

# net = network_main('/workspaces/PROKODE-DOCKER', '/workspaces/PROKODE-DOCKER/src/inputs/annotation.csv', '/workspaces/PROKODE-DOCKER/src/tfbs.csv', '/workspaces/PROKODE-DOCKER/src/decay_rates.csv')
