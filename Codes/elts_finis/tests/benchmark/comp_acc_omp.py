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
    print("Compiler:        " + compiler[i] + "\nMail size:     " + str(m_size[i]) + "\nTime elapsed:   " + str(time[i]))


def gflops(d):
    r = []
    for i in range(len(d["time"])):
        r.append(4 *  d["m_size"][i]/ (d["time"][i] * (10 ** 9)))
    return r

######### --- Raw times graph --- ########

plt.plot(perf["pgcc"]["m_size"], perf["pgcc"]["time"], color='blue', label="Parallel")
plt.plot(perf["icc"]["m_size"], perf["icc"]["time"], color='orange', label="OpenMP")
plt.legend()
plt.title("Comparaison des temps d'exécution d'implémentations du solveur de Jacobi")
plt.xlabel("Taille du maillage ($10^{-3}$ unités arbitraires)")
plt.ylabel("Temps d'exécution (secondes)")

######## --- GFLOPS graph --- ########

# plt.plot(perf["acc"]["m_size"], gflops(perf["acc"]), color='blue', label="Parallel")
# plt.plot(perf["ker"]["m_size"], gflops(perf["ker"]), '--', color='blue', label="Kernels")
# plt.plot(perf["omp"]["m_size"], gflops(perf["omp"]), color='orange', label="OpenMP")
# plt.legend()
# plt.title("Comparaison des performances en GFLOPS d'implémentations de l'intégration\n par la méthode des trapèzes")
# plt.xlabel("Nombre de trapèzes")
# plt.ylabel("Nombre de GFLOPS")

plt.show()

