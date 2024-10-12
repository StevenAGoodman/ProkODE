from multiprocessing import Pool
import linecache
import time
import json
from model_accuracy_functions import *

# inputs
TF_BETA_FILEPATH = "./tf_betas.json"

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
    return data_df, gene_key, tf_key, groups

def test_config_wBetas():
    while True:
      line = linecache.getline(TF_BETA_FILEPATH, latest_line)
      if line != "":
        break
      time.sleep(2)

    line = json.loads(line)
    tf_key =
    beta_values = list(line["beta values"].values())

    test_config_accuracy()



def call_py(feature_string, training_df, gene_key, tf_key,groups, prokode_dir):
    time.sleep(1)
    print(feature_string)
    print(training_df)
    network_loc = prokode_dir + "/src/network.json"
    with open(PYSCRIPT_FILE_PATH, 'r') as pyscript:
        data_file = DATA_FILE_PATH
        main(feature_string, training_df, prokode_dir, network_loc, gene_key, tf_key, groups, data_file)

if __name__ == '__main__':
    network_loc = PROKODE_BASEDIR + "/src/network.json"
    training_df, gene_key, tf_key, groups = init_stuffs(DATA_FILE_PATH, PROKODE_BASEDIR + "/src/network.json")
    training_df[]

    pool = Pool(processes=15)
    results = [pool.apply_async(call_py(x,training_df,gene_key, tf_key, groups,PROKODE_BASEDIR), (x,), callback=callback) for x in perms[1:10]]

    for result in results:
        result.wait()

