#generates the data matrix for random graphs after a threshold on probabilities
from __future__ import division
import scipy.integrate as spi
import numpy as np
import networkx as nx
import random
#import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set()

n = 1000
k = 0.1
pc = 0.6
pd = 0.3
pn = 0.1
#p = 0.0001

def graph_generate(p):
    G = nx.erdos_renyi_graph(n, p, seed=123, directed=False)
    return G

def graph_properties(G):    
    CC = nx.average_clustering(G, nodes=None, weight=None, count_zeros=True)
    spl = nx.average_shortest_path_length(G, weight=None)
    den = nx.density(G)
    dia = nx.diameter(G, e=None)
    deg = G.degree(nbunch=None, weight=None)
    avgdeg = np.mean(np.array(dict(deg).values()))
    max_deg = max(np.array(dict(deg).values()))
    return np.array([CC, spl, den, dia, avgdeg, max_deg])

prob = np.linspace(0.0072, 0.5, 5)                     
time = np.linspace(0, 100, 100)
Nodes = np.arange(0, n, 1)
Data_matrix = np.zeros((len(prob), 14))
count = -1
for p in range(len(prob)):
    print p
    S_init = random.sample(Nodes, 995)
    I_init = random.sample(np.delete(Nodes, S_init), 5)
    R_init = np.delete(Nodes, np.concatenate((S_init, I_init)))
    D = {}
    S = [len(S_init)]
    I = [len(I_init)]
    R = [len(R_init)]
    a = np.zeros(20)
    b = np.zeros(20)
    c = np.zeros(20)
    e = np.zeros(20)
    Ro = np.zeros(20)
    Ro_approx = np.zeros(20)
    Ro_approx2 = np.zeros(20)
    S_last_20 = np.zeros(20)
    I_last_20 = np.zeros(20)
    R_last_20 = np.zeros(20)
    G = graph_generate(prob[p])
    count = count + 1
    for t in range(len(time)):
        dc = nx.number_connected_components(G)
        if (dc > 1):
            print "break"
            break
        delta_I = 0
        delta_R = 0
        for i in range(len(G.nodes())):
            if (list(G.nodes())[i] in S_init):
                    D[list(G.nodes())[i]] = 0
            if (list(G.nodes())[i] in I_init):
                    D[list(G.nodes())[i]] = 1
            if (list(G.nodes())[i] in R_init):
                    D[list(G.nodes())[i]] = 2
        
        Suscept_count = len(S_init)
        Infect_count = len(I_init)
        Recov_count = len(R_init)
        
        for i in range(len(S_init)):
             neighbor = G.neighbors(S_init[i])
             Infect_neigh_count = 0
             for j in neighbor:
                 if (D[j] == 1):
                     Infect_neigh_count = Infect_neigh_count + 1
                     
             pv = (1-np.exp(-k*Infect_neigh_count))
             if (pv >= 0.09):
                 delta_I = delta_I + 1
                 Suscept_count = Suscept_count - 1
                 Infect_count = Infect_count + 1

        
        #if  (Infect_count >= 6):
        delta_R = int(pc * Infect_count)
        Recov_count = Recov_count + delta_R
        Infect_count = Infect_count - delta_R
        Suscept_count =  Suscept_count + int(pn* Recov_count)  + int(pd * Infect_count)
        Infect_count = Infect_count - int(pd * Infect_count)
        Recov_count = Recov_count - int(pn* Recov_count)
        #print Suscept_count, Infect_count, Recov_count          
        S_init = random.sample(Nodes, Suscept_count)
        I_init = random.sample(np.delete(Nodes, S_init), Infect_count)
        R_init = np.delete(Nodes, np.concatenate((S_init, I_init)))         
        S.append(Suscept_count)
        I.append(Infect_count)
        R.append(Recov_count)

        if (t > 79):
            A = delta_I/(Suscept_count*Infect_count)
            B = delta_R/ Infect_count
            C = (1- (delta_R/ Infect_count))*(int(pd * Infect_count)/Infect_count)
            a[t - 80] = A
            b[t - 80] = B
            c[t - 80] = C
            e[t - 80] = int(pn* Recov_count)/Recov_count
            Ro[t - 80] = (A*n)/(B+C)
            Ro_approx[t - 80] = A/B
            Ro_approx2[t - 80] = A/(B+C)
            S_last_20[t-80] = Suscept_count
            I_last_20[t-80] = Infect_count
            R_last_20[t-80] = Recov_count
                
                
    if (dc == 1):
        Ro_N = np.mean(Ro)
        Ro_approx_N = np.mean(Ro_approx)
        Ro_approx2_N = np.mean(Ro_approx2)
        S_N = np.mean(S_last_20)
        I_N = np.mean(I_last_20)
        R_N = np.mean(R_last_20)
        Ip = max(I)
        #IIp = int(np.where(np.array(I)==Ip)[0])-1
        a = np.array(I)
        IIp = np.unravel_index(np.argmax(a, axis = 0), a.shape)[0]
        P = graph_properties(G)
        Data_matrix[int(count)] = np.array([P[0], P[1], P[2], P[3], P[4], P[5], S_N, I_N, R_N, Ro_N, Ip, IIp, Ro_approx_N*P[4], Ro_approx2_N*P[4]])
print Data_matrix        
#np.save("Data_matrix_random.npy", Data_matrix)
