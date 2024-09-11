import numpy as np
import pandas as pd
import subprocess
import os

def run_CiiiDER(promoters_loc, matrices_loc, deficit_val, prokode_dir):
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
    subprocess.run(['java','-jar','C:/Users/cryst/OneDrive/Documents/LOFScreening/CiiiDER_TFMs/CiiiDER.jar', '-n', prokode_dir +'/src/preprocessing/config.ini'])
    c_output_fpath = c_output_fpath[:-3] + 'csv'

    # NEED TO GET RID OF LINE BREAKS
    return c_output_fpath

def create_promoterf(genome_loc, annotation_loc, prokode_dir):
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
    return np.exp(score)

def create_tfbs(ci_results_loc, prokode_dir):
    # import ciiider results
    ciiider_df = pd.read_csv(ci_results_loc, names=['tg','_','tf','tf_matrixid','start','end','strand','prescore','score','seq'])

    tfbs_df = ciiider_df[['tf','tg']].copy()

    kd_vals = []
    for score in ciiider_df['score'].to_list():
        kd = score_to_kd(score)
        kd_vals.append(kd)

    tfbs_df.insert(2, 'Kd', kd_vals)

    # write to tfbs.csv
    tfbs_df.to_csv(open(prokode_dir + '/src/tfbs.csv', 'w'))

def main(prokode_dir, genome_loc, annotation_loc, pfm_database_loc):
    # # create promoters.fa
    # promoters_loc = create_promoterf(genome_loc,annotation_loc, prokode_dir)

    # # run CiiiDER with files
    # CiiiDER_results_loc = run_CiiiDER(promoters_loc,pfm_database_loc, 0.1, prokode_dir)
    CiiiDER_results_loc = prokode_dir + '/src/preprocessing/CiiiDER_results.csv'
    # configer into tf binding site csv
    create_tfbs(CiiiDER_results_loc, prokode_dir)

# temp start code
main('C:/Users/cryst/OneDrive/Documents/LOFScreening/PROKODE','C:/Users/cryst/OneDrive/Documents/LOFScreening/PROKODE/src/inputs/genome.fasta','C:/Users/cryst/OneDrive/Documents/LOFScreening/PROKODE/src/inputs/annotation.csv','C:/Users/cryst/OneDrive/Documents/LOFScreening/PROKODE/src/preprocessing/pfmdb.txt')