import numpy as np
import pandas as pd
import random
import json

max_runtime = 100
# f = [] # r x 1
# A = [] # m x r
# P = [] # r x n
# E = [] # m x n

with open("P.json", "r") as file:
    P = np.matrix(np.array(json.loads(file.read())))
with open("f.json", "r") as file:
    f = np.matrix(np.array(json.loads(file.read())))
with open("A_init.json", "r") as file:
    A = np.matrix(np.array(json.loads(file.read())))



######## is there a python (AI?) package that can do this loss function?
######## Loss function: min|P - E * A^+|^2

#create A valid for least squares
A_inv = np.linalg.pinv(A)
# block_a = np.block([[A_inv,A_inv],[A_inv,A_inv]])
# print(block_a)

## Initail iteration
b = column_matrix = P.ravel(order='F').reshape(-1, 1)

E,_,_,_ = np.linalg.lstsq(A_inv,b) # solves min |b - ax|

check = E - (A @ b)
print(check)
# for t in range(max_runtime):

