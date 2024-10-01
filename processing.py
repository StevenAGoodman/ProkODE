import pandas as pd
import numpy as np
import json
import statistics
from src.run import create_network_json_main
import cvxpy 

# be sure to run run.py before running this file
network_loc =

# global jazz
num_data = 10
sensor_normal_dist = 10
basal_rate = 3
decay_rate = np.log(2)/300
Np = 6000
Kd_p = 0.1
Nns = 4600000

def get_data_unit(i):
    
    
# input: pair of time points

for i in range(num_data):
    data_unit = get_data_unit(i)

    
