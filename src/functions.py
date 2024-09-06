# grand software inputs: nucleotide sequence of entire genome: fasta, geneIDs & their start coordinates: csv (geneid,start)
import numpy as np
import pandas as pd
import json
import subprocess

downstream = 70
upstream = 30

def eq1(beta, N_tf, Kd_tf, N_p, Kd_p, Nns):
    P_tf = N_tf / (N_tf + Kd_tf)
    E_basal = N_p / (Nns * Kd_p)
    E_PTF = beta * E_basal
    R = P_tf * E_PTF + (1 - P_tf) * E_basal
    return R

def rev_eq1(R,N_tf, Kd_tf, N_p, Kd_p, Nns):
    P_tf = N_tf / (N_tf + Kd_tf)
    E_basal = N_p / (Nns * Kd_p)
    E_PTF = (R - ((1 - P_tf) * E_basal)) / P_tf
    beta = E_PTF / E_basal
    return beta

def get_pwm(tf):
    pwm_dataset = open('example/pwm_dataset.json')
    json.load(pwm_dataset)
    pwm_dataset.close()

def make_promoter_seq_file(genome_path,genome_an_path):
    df = pd.read_csv(genome_an_path, names=['gene', 'start'])
    genome_file = open(genome_path, 'r')
    genome = genome_file.read() 
    genome_file.close()

    fasta_txt = ''
    for row in len(df.index):
        gene = df.loc[row, 'gene']
        start = df.loc[row,'start']
        seq = genome[start - downstream: start + upstream]

        fasta_txt += f'>{gene}\n{seq}\n\n'

    file_name = 'example/promoters.fa'
    with open(file_name, 'w') as file:
        file.write(fasta_txt)
    return file_name

def get_kd(delta_g):
    return np.exp(delta_g / (R * T))

def run_trap(sequence_file, pwm_file, output_file):
    subprocess.run(['annotate', '-s', sequence_file, '-pwm', pwm_file, '-o', output_file])

def create_tfbs_annotations(tf_arr, score_thresh, genome_path, gene_annot_path):
    promo_seq_file = make_promoter_seq_file(genome_path,gene_annot_path)
    
    output = ''
    for tf in tf_arr:
        pwm_file = get_pwm(tf)
        run_trap(promo_seq_file, pwm_file, 'example/trap_output.txt')
        ## analyze trap for genes that get bind and affinites
        for gene in trapresults:
            output += f'{tf}\t{promoter_gene}\t{Kd}\n'

    with open('example/tfbs.tsv', 'w') as file:
        file.write('tf\tgene\tKd\n')
        file.write(output)

            


# all the data we need for a tf-rnap-tf interaction is the prob of rnap binding at instant, N_tf, Kd_tf, N_p, Kd_p, Nns
# how does gene expression relate to prob of rnap binding at instant??? 

# def subset_tf_of_genes(gene_arr):

#     return tf_arr


# def express_to_R(exp):
#     return R

# def get_tgs(tf):
#     return tg_arr

# def main():
#     df = pd.read_csv('data.csv', names=['gene', 'expression'])
#     tfs = subset_tf_of_genes(df['gene'].values)

#     Nns = 4600000
#     p_gene = ''
#     N_p = df.loc[p_gene, 'expression']

#     with open('results.csv', 'w') as file:
#         file.write('tf,tg,beta')
#         for tf in tfs:
#             N_tf = df.loc[tf, 'expression']
#             tg_arr = get_tgs(tf)
#             for tg in tg_arr:
#                 Kd_tf = get_Kd(tf, tg)
#                 Kd_p = get_Kd(p_gene, tg)
#                 R = express_to_R()
#                 beta = rev_eq1(R, N_tf, Kd_tf, N_p,Kd_p,Nns)
#                 file.write(f'{tf},{tg},{beta}')

