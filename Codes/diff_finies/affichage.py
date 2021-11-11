import csv
import matplotlib.pyplot as plt
import numpy as np
from math import *

iterations = []
size = []
execution = []


with open('results.csv', 'r') as csvfile:
    input_datas = csv.reader(csvfile, delimiter=',')
    for row in input_datas:
        iterations.append(int(row[1]))
        size.append(int(row[0]))
        execution.append(float(row[2]))

flops = [(iterations[i] * (size[i] ** 2) * 12) / (execution[i] * (10 ** 9)) for i in range(len(iterations))]

plt.plot(size, flops)
plt.title("Nombre de GFLOPS en fonction de la taille de la grille")
plt.xlabel("Taille de la grille")
plt.ylabel("GFLOPS")
plt.show()

