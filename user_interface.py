### UNDER DEVELOPMENT ###
from src.run import *

# paramters
reset = True
prokode_dir = '/workspaces/git'
data_file = 'GSE90743_E14R025_raw_counts.txt'
genome_loc = '/workspaces/git/src/inputs/genome.fasta'
annotation_loc = '/workspaces/git/src/inputs/annotation.csv'
pfm_database_loc = '/workspaces/git/src/preprocessing/pfmdb.txt'
CiiiDER_jar_loc = '/CiiiDER/CiiiDER_TFMs/CiiiDER.jar'
CiiiDER_thresh = 0.5

def user_interface_main():
    network_loc = create_network_json_main()
    create_interface_files(network_loc)