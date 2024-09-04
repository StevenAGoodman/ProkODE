import numpy as np
import pandas as pd

# set up functions
current_path = # user input
import sys
sys.path.append(f'{current_path}/src')
sys.path.append(f'{current_path}/src/preprocessing')
from config_promo import *
from get_tfbs_decays import *

config_promoters(genome_path, annotation_path)
## ...continue...
