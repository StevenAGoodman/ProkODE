from multiprocessing import Pool
import itertools
import subprocess
import random

PYTHON_EXEC_PATH = "/usr/local/bin/python"
PYSCRIPT_FILE_PATH = "/workspaces/ProkODE/beta_training/processing.py"
PROKODE_BASEDIR =  "/workspaces/ProkODE" # !! NO end '/' !!
DATA_FILE_PATH = "/workspaces/ProkODE/beta_training/GEO_expression_data/combined_gene_expression.csv"

chars_list = ["@AB", "@AB", "@A", "@AB", "@", "@", "@", "@A", "01", "01", "01", "01", "01", "01", "01", "01"]

chars_list = [list(chars) for chars in chars_list]

perms = []
for i in itertools.product(*chars_list):
    perms.append("".join(i))

random.shuffle(perms)

def call_py(feature_str):
    print(feature_str)
    subprocess.run([PYTHON_EXEC_PATH, PYSCRIPT_FILE_PATH, feature_str, DATA_FILE_PATH, PROKODE_BASEDIR])

p = Pool()
results = p.map(call_py, perms)

p.close()
p.join()