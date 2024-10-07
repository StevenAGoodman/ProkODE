import numpy as np
import pandas as pd
import re
import json

# constants... that may change
temperature = 298 # kelvin
elongation_rate = 60 # nt/s
peptide_rate = 20 # aa/s
Kd_ribo_mrna = 0.1 # WHAT IS THE BINDING AFFINITY OF RIBO TO DALGARNO SEQ???? maybe use molecular docking
len_taken_by_rnap = 30 # nt
len_taken_by_ribo = 30 # aa



# fundemental functions
def score_to_K(score):
    delta_G = -1 * score * 1000
    return 1 / np.exp(delta_G / (1.98722 * temperature))
def search_network_json(network_loc, gene_str:str):
    gene_str = str.lower(gene_str)
    with open(network_loc, 'r') as network_file:
        for _, line in enumerate(network_file):
            print(line)
            if len(line) > 3 and re.search(gene_str, str.lower(line[1:line.find('":{')])) != None:
                gene_info = line[line.find('":{')+2:-2].replace("\n",'')
                return json.loads(gene_info)
            elif len(line) > 3 and gene_str in json.loads(re.search('"synonyms": (\[.+\]), ', line).group(1)):
                return json.loads(line[line.find('":{')+2:-2].replace("\n",''))
            else:
                continue



#########            #########
### FUNCTIONS FOR TWEAKING ###
#########            #########



# max rates
R_max_txn = len_taken_by_rnap / elongation_rate
R_max_trans = len_taken_by_ribo / peptide_rate



# binding probability functions
def tf_probabiltiy(N_tf, Kd_tf_target, genome_len):
    return N_tf / (genome_len * (N_tf + Kd_tf_target))
def rnap_probabiltiy(N_rnap, Kd_rnap_target, genome_len):
    return N_rnap / (genome_len * (N_rnap + Kd_rnap_target))
def ribo_probabiltiy(N_ribo, ):
    return N_ribo / (N_ribo + Kd_ribo_mrna)



# Transcription rate function
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



# Beta function
def beta_function(P_arr, *beta_arr):
    assert len(P_arr) == len(beta_arr)
    return np.array(P_arr) @ np.log10(np.array(beta_arr)).T
def beta_from_context(R_txn, P_rnap):
    return R_txn / (P_rnap * R_max_txn)



# Translation rate function
def translation_rate(gene, protein_amnts, N_ribo, gene_info_dict):
    if gene_info_dict == None:
        return 0
    # calculate max translation rate
    mRNA_len = gene_info_dict["mRNA length"]
    max_translation_rate = 1 / (len_taken_by_ribo / peptide_rate)

    # calculate ribo binding prob
    P_ribo_bound = N_ribo / (N_ribo + Kd_ribo_mrna)

    return P_ribo_bound * max_translation_rate



# mRNA decay model
def RNA_decay_rate(gene, prev_total_mRNA_amnt, protein_amnts, decay_dict, gene_key):
    growth_rate = get_grow_rate(protein_amnts)

    half_life = 0

    # average half life from bionumbers
    half_life += (4.5148 + np.random.random()) * 60 # seconds

    # # calculation of natural decay component
    # natural_mRNA_half_life = 0

    # # get degrad_prots of gene
    # decay_info = decay_dict#[gene] # should be a dict of

    # rate_of_mRNA_cleavage = 0
    # for degrad_prot, K in decay_info.items():
    #     N_dp = protein_amnts[gene_key.index(degrad_prot)]
    #     # K_1 is the reaction rate from [Enzyme] + [mRNA] -> [Enzyme-mRNA] (ie, K_1[Enzyme][mRNA] = [Enzyme-mRNA] create / time ) ... k_1 is in 1 / (mols * time)
    #     rate_of_mRNA_cleavage += K * prev_total_mRNA_amnt * N_dp / (1 + K * prev_total_mRNA_amnt)
    #     # what is the relationship of multiple degrading proteins? it should be additive, right?
    #     # how can you identify a protien as degrading?
    # #     # what is the relationship of degrading prots component and decay (ie half life) component

    # half_life += rate_of_mRNA_cleavage + natural_mRNA_half_life

    return (np.log(2) / half_life) + growth_rate



# Protein decay model
def protein_decay_rate(gene, protein_amnts, decay_dict, gene_key):
    growth_rate = get_grow_rate(protein_amnts)

    half_life = 0 # mins
 
    # average half life from bionumbers
    half_life += (664.75 + np.random.random() * 100) * 60 # seconds

    # total_half_life = 0

    # # N end rule
    # total_half_life += 1 / decay_dict["N end"][gene]

    # # decay from misfolded proteins
    # total_half_life += 1 / decay_dict["misfold"][gene]()

    # return 1 / total_half_life
    return (np.log(2) / half_life) + growth_rate



# RNAP change model
def RNAP_amount(prev_rnap, protein_amnts):
    return 2200 # μm^3
    


# Ribosome change model
def ribo_amount(prev_ribo, protein_amnts):
    return 3000 # μm^3



# Sigma factor competition
def sigma_competition():
    return None



# Cell growth
def get_grow_rate(protein_amnts):
    growth_rate = 0
    return growth_rate



# link to processing.py
def beta_from_overall_mRNA(gene, gene_mRNA_amnt, overall_mRNA_change_rate, protein_amnts, regulators_dict:dict, gene_key, tf_key:list, genome_length, cell_volume, Kd_rnap, prev_rnap, prev_ribo):
    N_rnap = RNAP_amount(prev_rnap, protein_amnts)
    N_ribo = ribo_amount(prev_ribo, protein_amnts)

    mRNA_decay_rate = RNA_decay_rate(gene, gene_mRNA_amnt, protein_amnts, decay_dict, gene_key)
    mRNA_creation_rate = overall_mRNA_change_rate + mRNA_decay_rate * gene_mRNA_amnt

    P_rnap = rnap_probabiltiy(N_rnap, Kd_rnap, genome_length)
    beta_all = beta_from_context(mRNA_creation_rate, P_rnap)

    coefficient_arr = [0] * len(tf_key)
    for regulator, reg_details in regulators_dict.items():
        P_tf = tf_probabiltiy(protein_amnts[gene_key.index(regulator)], reg_details["Kd"], genome_length)
        coefficient_arr[tf_key.index(regulator)] = P_tf

    return coefficient_arr, beta_all, N_rnap, N_ribo