import numpy as np
import pandas as pd
from preprocessing.preprocessing import preprocessing_main
from config_network_structure import network_main
# from plot_system import plot_system
from decay.processing import decay_main
import os

# paramters
reset = True
prokode_dir = 'C:/Users/cryst/LOFScreening/archive/PROKODE'
data_file = 'GSE90743_E14R025_raw_counts.txt'
genome_loc = prokode_dir + '/src/inputs/genome.fasta'
annotation_loc = prokode_dir + '/src/inputs/annotation.tsv'
operons_loc = prokode_dir + '/src/inputs/operons.tsv'
pfm_database_loc = prokode_dir + '/src/preprocessing/pfmdb.txt'
CiiiDER_jar_loc = './CiiiDER_TFMs/CiiiDER.jar'
CiiiDER_thresh = 0.6

# global jazz
sensor_normal_dist = 10
basal_rate = 3
decay_rate = np.log(2)/300
Np = 6000
Kd_p = 0.1
Nns = 4600000

def create_network_json_main(prokode_dir, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, add_betas = False, reset=True):
    # reset file structure
    reset_files = ['./results.json','./run.log','./src/tfbs.csv','./src/network.json','./src/__pycache__/config_network_structure.cpython-312.pyc','./src/preprocessing/CiiiDER_results.bsl', './src/preprocessing/CiiiDER_results.csv','./src/preprocessing/CiiiDER_results.tsv','./src/preprocessing/config.ini','./src/preprocessing/motif_matrices.ml','./src/preprocessing/motif_matrices.txt','./src/preprocessing/promoters.fa','./src/preprocessing/__pycache__/preprocessing.cpython-312.pyc']
    for file in reset_files:
        try:
            os.remove(file)
        except:
            continue

    # create tfbs.csv
    print('\n\tpreprocessing inputs...')
    tfbs_loc =  preprocessing_main(prokode_dir, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, add_betas) # '/workspaces/PROKODE-DOCKER/src/tfbs.csv', '/workspaces/PROKODE-DOCKER/src/decay_rates.csv'
    # tfbs_loc = 'C:/Users/cryst/LOFScreening/archive/PROKODE/src/tfbs.csv'
    # create network.json
    print('\tconfiguring network structure file...')
    network_loc = network_main(prokode_dir, annotation_loc, operons_loc, tfbs_loc)

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

create_network_json_main(prokode_dir, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, False, reset)
