### UNDER DEVELOPMENT ###
import numpy as np
import pandas as pd
from Bio import motifs

def get_specificities(out_loc, motifs_loc):
    return None    

def decay_protein_classifier(prokode_dir, molecule, gene_key:list):
    # function: intake a list of genes and identify all genes whose proteins effect decay of the molecule, write their motifs to a file

    if molecule == "mRNA":
        motifs_loc = prokode_dir + "/src/decay/mRNA_motifs.txt"
        motifs_out_loc = prokode_dir + '/src/decay/mRNA_motifs_present.txt'
    
    elif molecule == "protein":
        motifs_loc = prokode_dir + "/src/decay/protein_motifs.txt"
        motifs_out_loc = prokode_dir + '/src/decay/protein_motifs_present.txt'

    else: 
        raise Exception('molecule not found')

    motif_arr = []
    record = motifs.parse(open(motifs_loc, 'r'), 'JASPAR')
    for m in record:
        # search gene_key for motif name
        is_present = True if m.name in gene_key else False

        if is_present:
            motif_arr.append(m)

    out = motifs.write(motif_arr, 'JASPAR')
    open(motifs_out_loc, 'w').write(out)
    
    return motifs_out_loc

def get_protdecay_json(gene_key):
    decay_json = {}

    # N end rule
    n_end_arr = {}
    for gene in gene_key:

        n_end_arr[gene]

    decay_json["N-end"] = n_end_arr

    # misfolding
    
    return decay_json

def decay_main(prokode_dir, gene_key):
    mRNA_motifs_loc = decay_protein_classifier(prokode_dir, "mRNA", gene_key)
    mRNA_decay_bs_loc = get_specificities(prokode_dir + '/src/decay/mRNA_decay.csv', mRNA_motifs_loc)

    protein_motifs_loc = decay_protein_classifier(prokode_dir, "protein", gene_key)
    protein_decay_bs_loc = get_specificities(prokode_dir + '/src/decay/protein_decay.csv', protein_motifs_loc)
    protein_decay_json = get_protdecay_json()
    
    return {}, {} # mRNA_decay_json, protein_decay_json