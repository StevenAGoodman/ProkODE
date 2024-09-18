import numpy as np
import pandas as pd
from src.preprocessing.preprocessing import preprocessing_main
from src.config_network_structure import network_main
import os

def reset_file_structure(prokode_dir):
    os.remove()

def create_network_json_main(prokode_dir, genome_loc, annotation_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, add_betas = False):
    # reset file structure
    os.remove('./results.json')
    os.remove('./run.log')
    os.remove('./__pycache__/run.cpython-312.pyc')
    os.rmdir('./__pycache__')
    os.remove('./src/decay_rates.csv')
    os.remove('./src/tfbs.csv')
    os.remove('./src/__pycache__/config_network_structure.cpython-312.pyc')
    os.rmdir('./src/__pycache__')
    os.remove('./src/preprocessing/CiiiDER_results.bsl')
    os.remove('./src/preprocessing/CiiiDER_results.csv')
    os.remove('./src/preprocessing/config.ini')
    os.remove('./src/preprocessing/motif_matrices.ml')
    os.remove('./src/preprocessing/motif_matrices.txt')
    os.remove('./src/preprocessing/promoters.fa')
    os.remove('./src/preprocessing/__pycache__/preprocessing.cpython-312.pyc')
    os.rmdir('./src/preprocessing/__pycache__')

    # create tfbs.csv and decay_rates.csv
    print('\tpreprocessing inputs...')
    tfbs_loc, decay_rates_loc =  preprocessing_main(prokode_dir, genome_loc, annotation_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, add_betas) # '/workspaces/PROKODE-DOCKER/src/tfbs.csv', '/workspaces/PROKODE-DOCKER/src/decay_rates.csv'

    # create network.json
    print('\tconfiguring network structure file...')
    network_loc = network_main(prokode_dir, annotation_loc, tfbs_loc, decay_rates_loc)

    return network_loc

def user_interface_main():
    #
    a = ''