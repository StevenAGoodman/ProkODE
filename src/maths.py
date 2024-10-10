import numpy as np
import pandas as pd
import linecache
import re
import json

# constants... that may change




# fundemental functions
def score_to_K(score, temperature):
    delta_G = -1 * score * 1000
    return 1 / np.exp(delta_G / (1.98722 * temperature))
def create_network_key(network_loc):
    with open(network_loc, "r") as netfile:
        network_key = []
        for line in netfile:
            if len(line) < 3:
                continue
            else:
                line_keys = [str.lower(line[1:line.find('":{')])]
                line_keys.extend(json.loads(re.search('"synonyms": (\\[.+\\]), ', line).group(1)))
                network_key.append(line_keys)

    return network_key
def search_network_json(network_loc, network_key, gene_str:str):
    try:
        gene_str = str.lower(gene_str)
        net_index = [i for i in range(len(network_key)) if gene_str in network_key[i]]
        line = linecache.getline(network_loc, net_index)
        line = line.replace("\n","")[line.find('":{')+2:-1]
        return json.loads(line)
    except:
        return None

    

    # with open(network_loc, 'r') as network_file:
    #     for _, line in enumerate(network_file):
    #         if len(line) > 3 and re.search(gene_str, str.lower(line[1:line.find('":{')])) != None:
    #             gene_info = line[line.find('":{')+2:-2].replace("\n",'')
    #             return json.loads(gene_info)
    #         elif len(line) > 3 and gene_str in json.loads(re.search('"synonyms": (\\[.+\\]), ', line).group(1)):
    #             return json.loads(line[line.find('":{')+2:-2].replace("\n",''))
    #         else:
    #             continue



#########            #########
### FUNCTIONS FOR TWEAKING ###
#########            #########



# max rates
def max_rates(len_taken_by_rnap, elongation_rate, len_taken_by_ribo, peptide_rate, gene_info_dict, feature_id):
    if feature_id == "@" or feature_id == "A":
        R_max_txn = len_taken_by_rnap / elongation_rate
        R_max_trans = len_taken_by_ribo / peptide_rate
    elif feature_id == "B":
        transcript_len = gene_info_dict["transcript length"]
        mRNA_len = gene_info_dict["mRNA length"]
        R_max_txn = elongation_rate / transcript_len
        R_max_trans = peptide_rate / mRNA_len
    else:
        TypeError(f"feature id {feature_id} not found for {max_rates.__name__}")

    return R_max_txn, R_max_trans



# binding probability functions
def tf_probabiltiy(N_tf, Kd_tf_target, genome_len, feature_id):
    if feature_id == "@":
        return N_tf / (genome_len * (N_tf + Kd_tf_target))
    elif feature_id == "A":
        return N_tf / (N_tf + Kd_tf_target)
    elif feature_id == "B":
        return None
    else:
        TypeError(f"feature id {feature_id} not found for {tf_probabiltiy.__name__}")

def rnap_probabiltiy(N_rnap, Kd_rnap_target, genome_len, feature_fid):
    if feature_fid == "@":
        return N_rnap / (genome_len * (N_rnap + Kd_rnap_target))
    elif feature_fid == "A":
        return N_rnap / (N_rnap + Kd_rnap_target)
    elif feature_fid == "B":
        return None
    


# Transcription rate function
def transcription_rate(gene, protein_amnts, gene_key, N_rnap, Kd_rnap, gene_info_dict, genome_len, feature_id):

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
    # max_txn_rate = 1 / (len_taken_by_rnap / elongation_rate) # transcripts/s

    # calculate basal rnap binding prob
    P_rnap_basal = N_rnap / (genome_len * (N_rnap + Kd_rnap))

    return beta_all * P_rnap_basal * max_txn_rate



# Beta function
def beta_all_from_context(R_txn, N_rnap, Kd_rnap, genome_len, R_max_txn, feature_id):
    if feature_id == "@":
        P_rnap = N_rnap / (genome_len * (N_rnap + Kd_rnap))
        return R_txn / (P_rnap * R_max_txn) 
    elif feature_id == "A":
        P_rnap = N_rnap / (N_rnap + Kd_rnap)
        return R_txn / (P_rnap * R_max_txn)
    elif feature_id == "B":
        return  R_txn / (N_rnap * R_max_txn / (Kd_rnap + 1))
    else:
        TypeError(f"feature id {feature_id} not found for {beta_all_from_context.__name__}")


def beta_curveFit_func(feature_id):

    if feature_id == "@":
        def beta_function(P_arr, *beta_arr):
            assert len(P_arr) == len(beta_arr)
            return np.array(P_arr) @ np.array(beta_arr).T
        
    elif feature_id == "A":
        def beta_function(P_arr, *beta_arr):
            assert len(P_arr) == len(beta_arr)
            return np.log10(np.array(P_arr)) @ np.log10(np.array(beta_arr)).T
        
    elif feature_id == "B":
        def beta_function(P_arr, *beta_arr):
            assert len(P_arr) == len(beta_arr)
            return np.array(P_arr) @ np.log10(np.array(beta_arr)).T
        
    elif feature_id == "C":
        def beta_function(P_arr, *beta_arr):
            assert len(P_arr) == len(beta_arr)
            return np.array(P_arr) @ np.log2(np.array(beta_arr)).T
    
    return beta_function


# Translation rate function
def translation_rate(gene, gene_mrna_amount, protein_amnts, N_ribo, gene_info_dict, R_max_trans, Kd_ribo_mrna, feature_id):
    
    if feature_id == "@":
        translation_rate_raw = (N_ribo / (N_ribo + Kd_ribo_mrna)) * R_max_trans
        return translation_rate_raw * gene_mrna_amount
    
    elif feature_id == "A":
        translation_rate_raw = (N_ribo / (N_ribo + Kd_ribo_mrna)) * N_ribo * R_max_trans
        return translation_rate_raw * gene_mrna_amount

    elif feature_id == "B":
        return (N_ribo * R_max_trans * gene_mrna_amount) / (Kd_ribo_mrna + gene_mrna_amount)
    else:
        TypeError(f"feature id {feature_id} not found for {translation_rate.__name__}")


# mRNA decay model
def RNA_decay_rate(prev_gene_mRNA, protein_amnts, gene_info_dict, gene_key, growth_fid, binary_feature_arr):
    decay_info = gene_info_dict["mRNA decay"]
    binary_feature_arr = list(map(int, binary_feature_arr))

    half_life = 0

    if binary_feature_arr[0] == 1:
        # average half life from bionumbers
        half_life += (4.5148 + np.random.random()) * 60 # seconds

    if binary_feature_arr[1] == 1:
        # calculation of natural decay component
        half_life += 0

    if binary_feature_arr[2] == 1:
        half_life += get_grow_rate(protein_amnts, growth_fid)

    if binary_feature_arr[3] == 1:
        rate_of_mRNA_cleavage = 0
        for degrad_prot, K in decay_info.items():
            N_dp = protein_amnts[gene_key.index(degrad_prot)]
            # K_1 is the reaction rate from [Enzyme] + [mRNA] -> [Enzyme-mRNA] (ie, K_1[Enzyme][mRNA] = [Enzyme-mRNA] create / time ) ... k_1 is in 1 / (mols * time)
            rate_of_mRNA_cleavage += K * prev_gene_mRNA * N_dp / (1 + K * prev_gene_mRNA)
            # how can you identify a protien as degrading?
        #     # what is the relationship of degrading prots component and decay (ie half life) component
        half_life += np.log(2) / rate_of_mRNA_cleavage
    
    if binary_feature_arr[4] == 1:
        # Michaelis–Menten 
        rate_of_mRNA_cleavage = 0
        for degrad_prot, K in decay_info.items():
            N_dp = protein_amnts[gene_key.index(degrad_prot)]
            # K_1 is the reaction rate from [Enzyme] + [mRNA] -> [Enzyme-mRNA] (ie, K_1[Enzyme][mRNA] = [Enzyme-mRNA] create / time ) ... k_1 is in 1 / (mols * time)
            rate_of_mRNA_cleavage += (() * prev_gene_mRNA) / (K + prev_gene_mRNA)
            # how can you identify a protien as degrading?
        #     # what is the relationship of degrading prots component and decay (ie half life) component
        half_life += np.log(2) / rate_of_mRNA_cleavage

    return np.log(2) / half_life



# Protein decay model
def protein_decay_rate(gene, protein_amnts, gene_info_dict, gene_key, growth_fid, binary_feature_arr):
    binary_feature_arr = list(map(int, binary_feature_arr))
    prev_gene_prot = protein_amnts[gene_key.index(gene)]
    decay_info = gene_info_dict["protein decay"]

    half_life = 0 # secs

    if binary_feature_arr[0] == 1:
        # average half life from bionumbers
        half_life += (664.75 + np.random.random() * 100) * 60 # seconds

    if binary_feature_arr[1] == 1:
        # calculation of natural decay component
        half_life += 0

    if binary_feature_arr[2] == 1:
        half_life += get_grow_rate(protein_amnts, growth_fid)

    # # N end rule
    if binary_feature_arr[3] == 1:
        half_life +=  gene_info_dict["N end rule"]

    if binary_feature_arr[3] == 1:
        rate_of_mRNA_cleavage = 0
        for degrad_prot, K in decay_info.items():
            N_dp = protein_amnts[gene_key.index(degrad_prot)]
            # K_1 is the reaction rate from [Enzyme] + [mRNA] -> [Enzyme-mRNA] (ie, K_1[Enzyme][mRNA] = [Enzyme-mRNA] create / time ) ... k_1 is in 1 / (mols * time)
            rate_of_mRNA_cleavage += K * prev_gene_prot * N_dp / (1 + K * prev_gene_prot)
            # how can you identify a protien as degrading?
        #     # what is the relationship of degrading prots component and decay (ie half life) component
        half_life += np.log(2) / rate_of_mRNA_cleavage
    
    if binary_feature_arr[4] == 1:
        # Michaelis–Menten 
        rate_of_mRNA_cleavage = 0

        for degrad_prot, K in decay_info.items():
            N_dp = protein_amnts[gene_key.index(degrad_prot)]
            # K_1 is the reaction rate from [Enzyme] + [mRNA] -> [Enzyme-mRNA] (ie, K_1[Enzyme][mRNA] = [Enzyme-mRNA] create / time ) ... k_1 is in 1 / (mols * time)
            rate_of_mRNA_cleavage += (() * prev_gene_prot) / (K + prev_gene_prot)

        half_life += np.log(2) / rate_of_mRNA_cleavage

    # # decay from misfolded proteins
    if binary_feature_arr[6] == 1:
        prob = 3 # prob of misfolding to occur
        rate = 10
        half_life += prob * rate 

    # return 1 / total_half_life
    return np.log(2) / half_life



# RNAP change model
def RNAP_amount(prev_rnap, protein_amnts, feature_id):
    if feature_id == "@":
        return 2200 # μm^3
    
    elif feature_id == "A":
        # chance of proteins coming together form rnap
        return None
    else:
        TypeError(f"feature id {feature_id} not found for {RNAP_amount.__name__}")



# Ribosome change model
def ribo_amount(prev_ribo, protein_amnts, feature_id):
    if feature_id == "@":
        return 3000 # μm^3
    
    elif feature_id == "A":
        # chance of proteins coming together form ribo
        return None
    else:
        TypeError(f"feature id {feature_id} not found for {ribo_amount.__name__}")



# Sigma factor competition
def sigma_competition():
    return None



# Cell growth
def get_grow_rate(protein_amnts, feature_id):
    if feature_id == "0":
        growth_rate = 0

    elif feature_id == "A":
        growth_rate = 1/3600

    elif feature_id == "B":
        # check for genes correlated with cell growth
        growth_rate = None
    else:
        TypeError(f"feature id {feature_id} not found for {get_grow_rate.__name__}")

    
    return growth_rate



# link to processing.py
def beta_from_overall_mRNA():

    mRNA_decay_rate = RNA_decay_rate()
    mRNA_creation_rate = overall_mRNA_change_rate + mRNA_decay_rate * gene_mRNA_amnt

    beta_all = beta_all_from_context()

    coefficient_arr = [0] * len(tf_key)
    for regulator, reg_details in gene_info_dict["regulators"].items():
        if regulator == "polymerase":
            continue
        try:
            P_tf = tf_probabiltiy()
        except:
            # print(f"\t\t\tFailed: {regulator}\tgenome tf is not found in data's genes")
            continue
        coefficient_arr[tf_key.index(regulator)] = P_tf

    return coefficient_arr, beta_all

