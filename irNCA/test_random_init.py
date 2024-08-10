import numpy as np
import random
import json

n_tfs = 100
n_genes = 1000

A = [] # each inner array represent the cs for each tf, size: genes x tfs 
f = [] # percent of each tf active (range: 0 to 1)
P_1 = []

## Define f
# with open("f.json", "r") as f_file:
#     f = json.loads(f_file.read())
for tf in range(n_tfs):
    f.append(round((random.random()),2))

## Define A
for gene in range(n_genes):
    AInnerArr = []
    for tf in range(n_tfs):
        AInnerArr.append(round((random.random() * 10))) # !!improve function of probabilities
    A.append(AInnerArr)

## Define P_1
for tf in range(n_tfs):
    P_1.append(random.randint(0, 10))

## init data files
with open("A_init.json", "w") as A_file:
    A_file.write(json.dumps(A))
with open("f.json", "w") as f_file:
    f_file.write(json.dumps(f))
with open("P.json", "w") as P_file:
    P_file.write(json.dumps(P_1))