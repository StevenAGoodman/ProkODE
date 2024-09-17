import numpy as np
import pandas as pd
from src.preprocessing.preprocessing import preprocessing_main
from src.config_network_structure import network_main


def create_network_json_main(prokode_dir, genome_loc, annotation_loc, pfm_database_loc, CiiiDER_jar_loc, add_betas = False):
    # create tfbs.csv and decay_rates.csv
    tfbs_loc, decay_rates_loc =  preprocessing_main(prokode_dir, genome_loc, annotation_loc, pfm_database_loc, CiiiDER_jar_loc, add_betas) # '/workspaces/PROKODE-DOCKER/src/tfbs.csv', '/workspaces/PROKODE-DOCKER/src/decay_rates.csv'

    # create network.json
    network_loc = network_main(prokode_dir, annotation_loc, tfbs_loc, decay_rates_loc)

    return network_loc

def user_interface_main():
    #
    a = ''