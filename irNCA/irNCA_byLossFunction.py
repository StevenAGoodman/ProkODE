import numpy as np
import scipy 
import random
import json

import scipy.optimize

max_runtime = 100
n = 1000 # number of genes 
m = 1 # time is 1 at init
r = 100 # number of tfs
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

## Loss function: min|P - A^-1 * E|^2

#create A valid for least squares
A_inv = np.linalg.pinv(A)
# block_a = np.block([[A_inv,A_inv],[A_inv,A_inv]])
# print(block_a)

## Initail iteration
P = column_matrix = P.ravel(order='F').reshape(-1, 1)

huggabaloo = m*n+r*m

### Loss function optimization
def errorFunc(e,a):
    return np.linalg.norm(P - (a @ e))**2

def toVector(e, a):
    assert e.shape == (m, n)
    assert a.shape == (r, m)
    return np.hstack([e.flatten(), a.flatten()])

def toWZ(vec):
    assert vec.shape == (huggabaloo,)
    return vec[:m*n].reshape(m,n), vec[m*n:].reshape(r,m)

def doOptimization(e0,a0):
    def f(x): 
        e, a = toWZ(x)
        return errorFunc(e, a)

    result = scipy.optimize.minimize(f, toVector(e0, a0))
    # Different optimize functions return their
    # vector result differently. In this case it's result.x:
    result.x = toWZ(result.x) 
    return result

print(doOptimization(np.ones((m,n)), np.ones((r,m))))




