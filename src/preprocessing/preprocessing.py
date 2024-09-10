import numpy as np
import pandas as pd
import subprocess

def run_CiiiDER(promoters_loc, matrices_loc, deficit_val):
    # CiiiDER (https://ciiider.erc.monash.edu/) is a software that searches the promoter DNA regions of each gene with the binding motifs of each transcription factor to determine their binding sites. 

    c_output_fpath = '/workspaces/PROKODE-DOCKER/src/preprocessing/CiiiDER_results'
    
    # config config.ini
    config_o = f"""[General]
STARTPOINT = 1
ENDPOINT = 1\n
[Scan]
GENELISTFILENAME = {promoters_loc}
MATRIXFILE = {matrices_loc}
GENESCANRESULTS = {c_output_fpath}
DEFICIT = {deficit_val}"""
    
    open('/workspaces/PROKODE-DOCKER/src/preprocessing/config.ini', 'w').write(config_o)

    # run ciiider
    subprocess.run(['java','-jar','/CiiiDER/CiiiDER_TFMs/CiiiDER.jar', '-n', '/workspaces/PROKODE-DOCKER/src/preprocessing/config.ini'])
    return c_output_fpath

def create_promoterf(genome_loc, annotation_loc):
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
    promoters_loc = '/workspaces/PROKODE-DOCKER/src/preprocessing/promoters.fa'
    open(promoters_loc, 'w').write(promoter_contents)
    return promoters_loc

def score_to_kd(score):
    return np.exp(score)

def create_tfbs(ci_results_loc):
    # import ciiider results
    ciiider_df = pd.read_csv(ci_results_loc, names=['tg','_','tf','tf_matrixid','start','end','strand','prescore','score','seq'])

    tfbs_df = ciiider_df[['tf','tg']].copy()

    kd_vals = []
    for score in ciiider_df['score'].to_list():
        kd = score_to_kd(score)
        kd_vals.append(kd)

    tfbs_df.insert(2,kd_vals)

    # write to tfbs.csv
    tfbs_df.to_csv(open('/workspaces/PROKODE-DOCKER/src/tfbs.csv', 'w'))

def main(genome_loc, annotation_loc, pfm_database_loc):
    # create promoters.fa
    promoters_loc = create_promoterf(genome_loc,annotation_loc)

    # run CiiiDER with files
    CiiiDER_results_loc = run_CiiiDER(promoters_loc,pfm_database_loc, 0.1)

    # configer into tf binding site csv
    create_tfbs(CiiiDER_results_loc)

# temp start code
main('/workspaces/PROKODE-DOCKER/src/inputs/genome.fasta','/workspaces/PROKODE-DOCKER/src/inputs/annotation.csv','/workspaces/PROKODE-DOCKER/src/preprocessing/pfmdb.txt')