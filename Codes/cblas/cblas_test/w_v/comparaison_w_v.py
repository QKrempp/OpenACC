import csv
import matplotlib.pyplot as plt
import numpy as np
from math import *
    
gangs = []
workers = []
vectors = []
matrix = []
time = []

with open('results.csv', 'r') as csvfile:
    input_datas = csv.reader(csvfile, delimiter=',')
    for row in input_datas:
        gangs.append(int(row[0]))
        workers.append(int(row[1]))
        vectors.append(int(row[2]))
        matrix.append(int(row[3]))
        time.append(float(row[4]))

perf = {}

for i in range(len(gangs)):
    key = str(gangs[i]) + "g_" + str(workers[i]) + "w_" + str(vectors[i]) + "v"
    if key in perf:
        perf[key]["matrix_size"].append(matrix[i])
        perf[key]["time"].append(time[i])
        perf[key]["log_m_size"].append(log10(matrix[i]))
        perf[key]["log_time"].append(log10(time[i]))
    else:
        perf[key] = {"matrix_size":[matrix[i]], "time":[time[i]], "log_m_size":[log10(matrix[i])], "log_time":[log10(time[i])], "gangs":gangs[i], "workers":workers[i], "vectors":vectors[i]}

# ######### --- Log graph (Full) --- ########
# 
# plt.subplot(2, 2, 3)
# for i in perf.keys():
#     plt.plot(perf[i]["log_m_size"], perf[i]["log_time"], label=i)
# plt.legend()


# ######## --- Log Graph (Best) --- ########
# 
# plt.subplot(2, 2, 3)
# log_m_size = perf["512g_0w_0v"]["log_m_size"]
# z1 = np.polyfit(perf["512g_0w_0v"]["log_m_size"], perf["512g_0w_0v"]["log_time"], 1)
# p1 = np.poly1d(z1)
# z2 = np.polyfit(perf["512g_32w_32v"]["log_m_size"], perf["512g_32w_32v"]["log_time"], 1)
# p2 = np.poly1d(z2)
# plt.plot(perf["512g_0w_0v"]["log_m_size"], perf["512g_0w_0v"]["log_time"], 'o')
# plt.plot(log_m_size, p1(log_m_size))
# plt.plot(perf["512g_32w_32v"]["log_m_size"], perf["512g_32w_32v"]["log_time"], 'o')
# plt.plot(log_m_size, p2(log_m_size), label="modèle logarithmique")


######## --- Moyennes des rapports des temps d'axécution --- ########

plt.subplot(1, 1, 1)
hist = []
names = []
l = len(perf["512g_0w_0v"]["time"])
for i in perf.keys():
    if i != "512g_0w_0v":
        hist.append((sum([perf["512g_0w_0v"]["time"][j]/perf[i]["time"][j] for j in range(l)])/l) - 1)
        names.append("num_workers(" + str(perf[i]["workers"]) + ")\n vector_length(" + str(perf[i]["vectors"]) + ")")
x = np.arange(len(names))
names = [n for _,n in sorted(zip(hist, names))]
hist = sorted(hist)
plt.bar(x, hist, bottom=1)
plt.xticks(x, names, rotation=30, ha="center", fontsize=9)
plt.title("Moyenne des rapports de temps d'exécution\n avec le temps de référence (en secondes)")


# ######### --- Coefficient multiplicateur du temps d'exécution --- ########
# 
# plt.subplot(2, 2, 4)
# l = len(perf["512g_0w_0v"]["log_time"])
# for i in perf.keys():
#     if i != "512g_0w_0v":
#         sub=[10 ** (perf[i]["log_time"][j] - perf["512g_0w_0v"]["log_time"][j]) for j in range(l)]
#         plt.plot(perf[i]["log_m_size"], sub)


######## --- Modélisation --- ########
plt.show()
