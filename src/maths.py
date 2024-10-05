import numpy as np
import pandas as pd

# constants
temperature = 298 # kelvin
elongation_rate = 60 # nt/s
peptide_rate = 20 # aa/s
Kd_ribo_mrna = 0.1 # WHAT IS THE BINDING AFFINITY OF RIBO TO DALGARNO SEQ????
len_taken_by_rnap = 30 # nt
len_taken_by_ribo = 30 # aa

def score_to_K(score):
    delta_G = -1 * score * 1000
    return 1 / np.exp(delta_G / (1.98722 * temperature))

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

def creation_rate_from_overall_mRNA_stats(overall_mRNA_change_rate, mRNA_amnts, protein_amnts, gene_key):
    mRNA_decay_rate = RNA_decay_rate(gene, mRNA_amnts, protein_amnts, decay_dict, gene_key)

    return
