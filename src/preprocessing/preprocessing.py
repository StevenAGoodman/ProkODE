import numpy as np
import pandas as pd
import subprocess
import os
import sys

def run_CiiiDER(prokode_dir, promoters_loc, matrices_loc, deficit_val):
    # CiiiDER (https://ciiider.erc.monash.edu/) is a software that searches the promoter DNA regions of each gene with the binding motifs of each transcription factor to determine their binding sites. 

    c_output_fpath = prokode_dir + '/src/preprocessing/CiiiDER_results.txt'
    
    # config config.ini
    config_o = f"""[General]
STARTPOINT = 1
ENDPOINT = 1\n
[Scan]
GENELISTFILENAME = {promoters_loc}
MATRIXFILE = {matrices_loc}
GENESCANRESULTS = {c_output_fpath}
DEFICIT = {deficit_val}"""
    
    open(prokode_dir + '/src/preprocessing/config.ini', 'w').write(config_o)

    # run ciiider
    subprocess.run(['java','-jar','/CiiiDER/CiiiDER_TFMs/CiiiDER.jar', '-n', prokode_dir + '/src/preprocessing/config.ini'])

    c_output_fpath = c_output_fpath[:-3] + 'csv'

    return c_output_fpath

def create_promoterf(prokode_dir, genome_loc, annotation_loc):
    # config files
    genome = open(genome_loc, 'r').read()
    annotation_df = pd.read_csv(annotation_loc)

    # write promoter region surrounding each gene's start loc to promoters.fa
    promoter_contents = ''
    for i in annotation_df.index:
        gene = annotation_df.loc[i, 'geneid']
        start = annotation_df.loc[i, 'start']
        promo_seq = genome[start-150:start+50]
        promoter_contents += f'>{gene}\n{promo_seq}\n'

    # write to promoters.fa
    promoters_loc = prokode_dir + '/src/preprocessing/promoters.fa'
    open(promoters_loc, 'w').write(promoter_contents)
    return promoters_loc

def score_to_kd(score):
    return np.exp(-1 * score)

def add_beta_col(tfbs_loc):
    tfbs_df = pd.read_csv(tfbs_loc)

    beta_col = [np.random.randint(0,100) for i in tfbs_df.index]

    tfbs_df.insert(3,'beta', beta_col)

def create_tfbs(prokode_dir, ci_results_loc):
    # import ciiider results
    ciiider_df = pd.read_csv(ci_results_loc, names=['tg','_','tf','tf_matrixid','start','end','strand','prescore','score','seq'])

    tfbs_df = ciiider_df[['tf','tg']].copy()

    kd_vals = []
    for score in ciiider_df['score'].to_list():
        kd = score_to_kd(score)
        kd_vals.append(kd)

    tfbs_df.insert(2,'Kd', kd_vals)

    # write to tfbs.csv
    tfbs_loc = prokode_dir + '/src/tfbs.csv'
    tfbs_df.to_csv(open(tfbs_loc, 'w'))
    
    return tfbs_loc

def get_decay_rates(gene):
    return np.log(2)/300

def decay_rates_main(prokode_dir, annotation_loc):
    # read number of genes
    annotation_df = pd.read_csv(open(annotation_loc, 'r'))

    out = 'gene,decay_rate\n'
    for gene in annotation_df['geneid'].to_list():
        decay_rate = get_decay_rates(gene)
        out += f'{gene},{decay_rate}\n'

    decay_rates_loc = prokode_dir + '/src/decay_rates.csv'
    open(decay_rates_loc, 'w').write(out)

    return decay_rates_loc

def preprocessing_main(prokode_dir, genome_loc, annotation_loc, pfm_database_loc, add_betas=False):
    # create promoters.fa
    promoters_loc = create_promoterf(prokode_dir, genome_loc,annotation_loc)

    # run CiiiDER with files
    CiiiDER_results_loc = run_CiiiDER(prokode_dir, promoters_loc,pfm_database_loc, 0.1)

    # configer into tf binding site csv
    tfbs_loc = create_tfbs(prokode_dir, CiiiDER_results_loc)

    # config add betas if not trying to extract them
    if add_betas == True:
        add_beta_col(tfbs_loc)

    # create decay_rates.csv
    decay_rates_loc = decay_rates_main(prokode_dir, annotation_loc)
    
    return tfbs_loc, decay_rates_loc

# temp start code
preprocessing_main('/workspaces/PROKODE-DOCKER','/workspaces/PROKODE-DOCKER/src/inputs/genome.fasta','/workspaces/PROKODE-DOCKER/src/inputs/annotation.csv','/workspaces/PROKODE-DOCKER/src/preprocessing/pfmdb.txt')