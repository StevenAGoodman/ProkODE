from multiprocessing import Pool
import itertools
import subprocess
import random

PYTHON_EXEC_PATH = "/onfs/homes/01/a7hp6zz/python/3.10.10-venv/bin/python"
PROKODE_BASEDIR =  "/home/a7hp6zz/Repos/ProkODE" # !! NO end '/' !!
PYSCRIPT_FILE_PATH = f"{PROKODE_BASEDIR}/beta_training/processing.py"
DATA_FILE_PATH = f"{PROKODE_BASEDIR}/beta_training/GEO_expression_data/combined_gene_expression.csv"
num_procs = 32

chars_list = ["@AB", "@AB", "@A", "@AB", "@", "@", "@", "@A", "01", "01", "01", "01", "01", "01", "01", "01"]

chars_list = [list(chars) for chars in chars_list]

perms = []
for i in itertools.product(*chars_list):
    perms.append("".join(i))

random.shuffle(perms)

def call_py(feature_str):
    print(feature_str)
    subprocess.run([PYTHON_EXEC_PATH, PYSCRIPT_FILE_PATH, feature_str, DATA_FILE_PATH, PROKODE_BASEDIR])

p = Pool(num_procs)
results = p.map(call_py, perms)

p.close()
p.join()
