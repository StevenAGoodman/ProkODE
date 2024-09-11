import numpy as np
import pandas as pd
from src.preprocessing.preprocessing import *
from src.config_network_structure import *
from src.config_interface import *

def create_network_json_main(prokode_dir, genome_loc, annotation_loc):
    # create tfbs.csv and decay_rates.csv
    tfbs_loc, decay_rates_loc = preprocessing_main(prokode_dir, genome_loc, annotation_loc)

    # create network.json
    network_loc = network_main(prokode_dir, tfbs_loc, decay_rates_loc)

    return network_loc

def user_interface_main():
    #
    a = ''