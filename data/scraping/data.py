import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

def update(P_predict, H, R, z_k, x_predict):
    # kalman gain
    
    K_k = P_predict @ H.transpose() @ np.linalg.inv(H @ P_predict @ H.transpose() + R)
    
    # update
    x_est = x_predict + K_k @ (z_k - H @ x_predict)
    P_k = P_predict - K_k @ H @ P_predict

    return x_est, P_k

def predict(x_prevk, P_prevk, A, Q):
    x_predict = A @ x_prevk
    P_predict = (A @ P_prevk @ A.transpose()) + Q
    return x_predict, P_predict

lines = 1
max_points = -1

df = pd.read_csv("GSE90743_E14R025_raw_counts.txt",delimiter='\t', usecols=range(12))

# ploooooot stuff
t = [5,10,25,45,75,120,210,330,1500,1560,1680][:max_points] # time points
def func(t): return np.log(t) * 5 + np.random.uniform(-1.,1.) * 4

n = 2
plot1 = plt.subplot2grid((3, 3), (0, 0), colspan=n, rowspan=n)
plot2 = plt.subplot2grid((3, 3), (2, 0), rowspan=n, colspan=n)

for i in range(lines):
    y_i = df.iloc[i].values.flatten().tolist()[1:max_points]
    # y_i = [func(i) for i in t]
    plot1.plot(t, y_i, random.choice([ '-y', '-y']))

# kalman
x = np.array([[y_i[0]],
              [1.]]) # init state
prev_x = 0
P = np.array([[100.,0.],
              [0., 100.]]) # init uncertainty
A = np.array([[1.,.2],
              [0.,1.]])
H =  np.array([[1.,0.]])
Q = 0.01
R = 16

position = []
velocity = []

for k in range(len(t)):
    z = y_i[k]
    x, P = predict(x, P, A, Q)
    x, P = update(P, H, R, z, x)

    dx = (t[k] - prev_x)
    A = np.array([[1.,dx],
                  [0.,1.]])
    prev_x = t[k]

    position.append(x[0])
    velocity.append(x[1])


plot1.plot(t, position, '-r')
plot2.plot(t,velocity, '-g' )

plt.tight_layout()
plt.show()