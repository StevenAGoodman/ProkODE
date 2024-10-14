NUM_CORES = 6
PYTHON_EXEC_PATH = "python"
TF_BETA_FILEPATH = "c:/Users/cryst/LOFScreening/archive/PROKODE/ProkODE/partially_run_results.json"
DATA_FILE_PATH = "c:/Users/cryst/LOFScreening/archive/PROKODE/ProkODE/beta_training/GEO_expression_data/combined_gene_expression.csv"
PROKODE_BASEDIR =  "c:/Users/cryst/LOFScreening/archive/PROKODE/ProkODE" # !! NO end '/' !!

# libs
import multiprocessing as mp
import time
import pandas as pd
import numpy as np
from beta_estimation_functions import *
from model_accuracy_functions import *


latest_line = 1

def init_stuffs(data_file, network_loc):
    print("reading datafile to dataframe")
    # get data
    data_df = pd.read_csv(f"{data_file}", index_col="Unnamed: 0").multiply(1 / 6.02e23)
    # normalize as concentations
    # data_df = data_df.multiply(1 / 6.02e23)
    # define keys
    gene_key = list(data_df.index) # may have to subtracted header
    tf_key = [ tf for val in json.load(open(network_loc, 'r')).values() for tf in val["regulators"].keys() ]
    tf_key = list(set(tf_key)) # filter repeats
    return data_df, gene_key, tf_key

def getGroups(data_df):
    # filter out groups
    print("\t\t\t\tfiltering different groups")
    groups = [ [] for i in range(1000) ]
    for i in range(len(data_df.axes[1])):
        colname = data_df.columns[i]
        group_n = re.search("\\| (\\d+)", colname)
        if group_n == None:
            group_n = "0"
        else:
            data_df = data_df.rename(columns = {colname: colname[:colname.find(" | ")]})
            group_n = group_n.group(1)
        groups[int(group_n)].append(i)
    groups = [x for x in groups if x != []] # filter empties
    return data_df, groups




network_loc = PROKODE_BASEDIR + "/src/network.json"
training_df, gene_key, tf_key = init_stuffs(DATA_FILE_PATH, PROKODE_BASEDIR + "/src/network.json")
training_df = training_df.iloc[:, -20:]
training_df, groups = getGroups(training_df)
max_lines = 2594
network_key = create_network_key(network_loc)


done = []

def main(iter):
    global latest_line
    # if iter % 10 == 0:
    print("line:", iter)

    line = linecache.getline(TF_BETA_FILEPATH, iter)

    line = json.loads(line)
    config_id = [*line.keys()][0]
    tf_key = list([*line.values()][0]["beta values"].keys())
    beta_values = list([*line.values()][0]["beta values"].values())

    if type(beta_values[0]) == str:
     open("./test_completed.txt", "a").write(str(iter) + "\n")
     return "NA"     

    actual_matrix, predicted_matrix, gene_key, time_points_key = test_config_accuracy(training_df, tf_key, beta_values, network_loc, network_key, config_id)

    # do accuracy measurements
    general_acc = accuracy_representation(config_id, actual_matrix, predicted_matrix, gene_key, time_points_key, PROKODE_BASEDIR, True)
    if general_acc == np.float64("nan"):
        general_acc = str(general_acc) + ": floating point likely to small"
    open("./beta_training/config_results/overall.csv", 'a').write(f"{config_id},{general_acc}\n")

    open("./test_completed.txt", "a").write(str(iter) + "\n")


main(84)
# if __name__ == '__main__':
#     pool = mp.Pool(processes=NUM_CORES)
#     print("cpu count:", mp.cpu_count())
#     results = pool.map(main, list(range(1,max_lines)))
#     pool.close()
#     pool.join()
