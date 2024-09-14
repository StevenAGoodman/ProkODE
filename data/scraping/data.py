import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

df = pd.read_csv("GSE90743_E14R025_raw_counts.txt",delimiter='\t', usecols=range(12))
# 5m,n10, 25, 45, 75, 120, 210, 330, 25h, 26, 28

# ploooooot stuff
x = [5,10,25,45,75,120,210,330,1500,1560,1680]

fig, ax = plt.subplots()

for i in range(5):
    y_i = df.iloc[i].values.flatten().tolist()[1:]
    
    ax.plot(x, y_i, random.choice(['-r', '-g', '-b', '-y']))

# ax.scatter(x, y2, label='Dataset 2',marker='^')
# ax.scatter(x, y3, label='Dataset 3',marker='s')
# ax.scatter(x, y4, label='Dataset 4')
plt.show()