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
    else:
        perf[key] = {"mean_m_size":[], "mean_time":[], "matrix_size":[matrix[i]], "time":[time[i]], "log_m_size":[], "log_time":[], "gangs":gangs[i], "workers":workers[i], "vectors":vectors[i]}

for i in perf.keys():
    for j in range(len(perf[i]["matrix_size"])):
        if perf[i]["matrix_size"][j] in perf[i]["mean_m_size"]:
            n = perf[i]["mean_m_size"].index(perf[i]["matrix_size"][j])
            perf[i]["mean_time"][n] += perf[i]["time"][j]
        else:
            perf[i]["mean_m_size"].append(perf[i]["matrix_size"][j])
            perf[i]["mean_time"].append(perf[i]["time"][j])

    n = len(perf[i]["matrix_size"])/len(perf[i]["mean_m_size"])
    for j in range(len(perf[i]["mean_m_size"])):
        perf[i]["mean_time"][j] /= n

for key in perf.keys():
    for i in perf[key]["mean_time"]:
        perf[key]["log_time"].append(log2(i))
        if key == "0g_32w_32v":
            print("temps: log(" + str(i) + ") = " + str(log2(i)))
    for i in perf[key]["mean_m_size"]:
        perf[key]["log_m_size"].append(log2(i))
        if key == "0g_32w_32v":
            print("taille: log(" + str(i) + ") = " + str(log2(i)))


# ######### --- Log graph (Full) --- ########

plt.subplot(1, 1, 1)
for i in perf.keys():
    if perf[i]["gangs"] < 128:
        plt.plot(perf[i]["log_m_size"], perf[i]["log_time"], label=i)
    elif perf[i]["gangs"] == 0:
        plt.plot(perf[i]["log_m_size"], perf[i]["log_time"], '-o', label=i)

plt.title("Temps d'exécution en fonction de la taille de la matrice\n (échelle logarithmique en base 2)")
plt.xlabel("Taille des matrices carrées multipliées (log2)")
plt.ylabel("Temps d'exécution (log2)")
plt.legend(prop={"size":10})


######### --- Log graph (Gflops) --- ########

# fig, ax1 = plt.subplots()
# 
# color = 'tab:red'
# 
# ax1.plot(perf["0g_32w_32v"]["log_m_size"], perf["0g_32w_32v"]["log_time"], '-o', label="Courbe de référence", color=color)
# # ax1.title("Temps d'exécution en fonction de la taille de la matrice\n (échelle logarithmique en base 2)")
# ax1.set_xlabel("Taille des matrices carrées multipliées (log2)")
# ax1.set_ylabel("Temps d'exécution (log2)", color=color)
# ax1.legend(prop={"size":10})
# 
# gflops = []
# for i in range(len(perf["0g_32w_32v"]["mean_m_size"])):
#     gflops.append((perf["0g_32w_32v"]["mean_m_size"][i] ** 2) * (3 * perf["0g_32w_32v"]["mean_m_size"][i] + 2) / ((perf["0g_32w_32v"]["mean_time"][i]) * (10 ** 9)))
# 
# ax2 = ax1.twinx()
# color = 'tab:blue'
# ax2.plot(perf["0g_32w_32v"]["log_m_size"], gflops, label="GFLOPS théoriques", color=color)
# ax2.set_ylabel("GFLOPS", color=color)
# ax2.legend(prop={"size":10})
# 
# plt.title("Temps d'exécution en fonction de la taille des matrices comparé aux GFLOPS")
# fig.tight_layout()


# plt.subplot(2, 2, 2)
# for i in perf.keys():
#     if perf[i]["gangs"] >= 128:
#         plt.plot(perf[i]["mean_m_size"], perf[i]["mean_time"], label=i)
#     elif perf[i]["gangs"] == 0:
#         plt.plot(perf[i]["mean_m_size"], perf[i]["mean_time"], '--', label=i)
# 
# plt.legend(prop={"size":6})
# 

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

# plt.subplot(2, 1, 2)
# hist = []
# names = []
# l = len(perf["0g_32w_32v"]["time"])
# for i in perf.keys():
#     if i != "0g_32w_32v" and perf[i]["gangs"] >= 128:
#         hist.append((sum([perf["0g_32w_32v"]["time"][j]/perf[i]["time"][j] for j in range(l)])/l) - 1)
#         names.append("num_gangs(" + str(perf[i]["gangs"]) + ")")
# x = np.arange(len(names))
# names = [n for _,n in sorted(zip(hist, names))]
# hist = sorted(hist)
# plt.bar(x, hist, bottom=1)
# plt.xticks(x, names, rotation=30, ha="center", fontsize=9)
# plt.title("Moyenne des rapports de temps d'exécution\n avec le temps de référence (en secondes)")


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
