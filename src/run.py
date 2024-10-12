import os
import sys
import glob
import numpy as np
import pandas as pd

# python script dependencies
from preprocessing.preprocessing import preprocessing_main
from config_network_structure import network_main
from plot_system import plot_system
from decay.processing import decay_main

# string to bool function
def to_bool(bool_str):
    # from https://stackoverflow.com/a/13562496/127837
    """Parse the string and return the boolean value encoded or raise an exception"""
    # replace basestring by str in python 3.x
    if isinstance(bool_str, str) and bool_str:
        if bool_str.lower() in ['true', 't', '1']: return True
        elif bool_str.lower() in ['false', 'f', '0']: return False

    #if here we couldn't parse it
    raise ValueError("%s is not recognized as a boolean value" % bool_str)

# paramters
prokode_dir = sys.argv[1].replace('\\', '/')
reset = to_bool(sys.argv[2])
print(prokode_dir, reset)

genome_loc = prokode_dir + '/src/inputs/genome.fasta'
annotation_loc = prokode_dir + '/src/inputs/annotation.tsv'
operons_loc = prokode_dir + '/src/inputs/operons.tsv'
pfm_database_loc = prokode_dir + '/src/preprocessing/pfmdb.txt'
# CiiiDER_jar_loc = './CiiiDER_TFMs/CiiiDER.jar'
CiiiDER_thresh = 0


def create_network_json_main(prokode_dir, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_thresh, add_betas = False, reset=True):
    print("==========================\nProkODE: Inferring Gene Regulatory Network (GRN)\n==========================\n")
    # reset file structure
    if reset:
        reset_files = ['/results.json','/run.log','/src/tfbs.csv','/src/preprocessing/CiiiDER_results.bsl', '/src/preprocessing/CiiiDER_results.csv','/src/preprocessing/CiiiDER_results.tsv','/src/preprocessing/config.ini','/src/preprocessing/motif_matrices.ml','/src/preprocessing/motif_matrices.txt','/src/preprocessing/promoters.fa']
        reset_files.extend([i.replace("\\", "/") for i in glob.glob(prokode_dir + "/**/*.pyc", recursive=True)])

        for file in reset_files:
            try:
                os.remove(f"{prokode_dir}/{file}")
                print("\tfile reset:", file)
            except:
                continue

    # create tfbs.csv
    print('\n\tpreprocessing inputs...')
    tfbs_loc, floating_genes =  preprocessing_main(prokode_dir, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_thresh, add_betas) # '/workspaces/PROKODE-DOCKER/src/tfbs.csv', '/workspaces/PROKODE-DOCKER/src/decay_rates.csv'
    # tfbs_loc = 'C:/Users/cryst/LOFScreening/archive/PROKODE/src/tfbs.csv'
    # create network.json
    print('\tconfiguring network structure file...')
    network_loc = network_main(prokode_dir, genome_loc, annotation_loc, operons_loc, tfbs_loc, floating_genes)

    # plot_system(prokode_dir, prokode_dir + '/src/network.json', 1, 1, 0.1,100)

    return network_loc

def create_interface_files(network_loc):
    ### UNDER DEVELOPMENT ###


    # genome_df = pd.read_csv('./inputs/annotation.csv')
    # gene_arr = list(set(genome_df['geneid'].to_list()))
    # # config_network_json(gene_arr,'decay_rates.csv', 'tfbs.csv')
    # params = json.load(open('network.json', 'r'))
    # s0 = [random.randint(0,100) for i in range(len(gene_arr))]
    # E_basal = 0.00434782608696
    # print(len(params))
    # plot_system(params, s0, gene_arr)

    return None

create_network_json_main(prokode_dir, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_thresh, False, reset)
