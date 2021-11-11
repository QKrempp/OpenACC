import csv
import matplotlib.pyplot as plt
import numpy as np
from math import *
    
compiler = []
m_size = []
time = []

with open('results.csv', 'r') as csvfile:
    input_datas = csv.reader(csvfile, delimiter=',')
    for row in input_datas:
        compiler.append(row[0])
        m_size.append(int(row[1]))
        time.append(float(row[2]))

perf = {}

for i in range(len(compiler)):
    if compiler[i] in perf:
        perf[compiler[i]]["m_size"].append(m_size[i])
        perf[compiler[i]]["time"].append(time[i])
    else:
        perf[compiler[i]] = {"m_size":[m_size[i]], "time":[time[i]]}
    print("Compiler:        " + compiler[i] + "\nMatrix size:     " + str(m_size[i]) + "\nTime elapsed:   " + str(time[i]))


def gflops(d):
    r = []
    for i in range(len(d["time"])):
        r.append((d["m_size"][i] ** 2) * (2 * d["m_size"][i] + 3) / (d["time"][i] * (10 ** 9)))
    return r

######### --- Raw times graph --- ########

# plt.plot(perf["acc"]["m_size"], perf["acc"]["time"], color='blue', label="Parallel")
# plt.plot(perf["ker"]["m_size"], perf["ker"]["time"], '--', color='blue', label="Kernels")
# plt.plot(perf["omp"]["m_size"], perf["omp"]["time"], color='orange', label="OpenMP")
# plt.legend()
# plt.title("Comparaison des temps d'exécution d'implémentations de DGEMM")

######## --- GFLOPS graph --- ########

plt.plot(perf["acc"]["m_size"], gflops(perf["acc"]), color='blue', label="Parallel")
plt.plot(perf["ker"]["m_size"], gflops(perf["ker"]), '--', color='blue', label="Kernels")
plt.plot(perf["omp"]["m_size"], gflops(perf["omp"]), color='orange', label="OpenMP")
plt.legend()
plt.title("Comparaison des performances en GFLOPS d'implémentations de DGEMM")

plt.show()

