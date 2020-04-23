#  This codes simulates SIR network on the ER network in mean field
# approximation; also it calculates R0 at each step and modifies model
# parameters accordingly; This starts recovery from first step itself.
# Also there is no threshold on the probability pv, beyond which infection of
# susceptible nodes starts.
from __future__ import division
import scipy.integrate as spi
import numpy as np
import networkx as nx
import random
import matplotlib.pyplot as plt
import seaborn as sns
import time
sns.set()
sns.set_style("white")
sns.set_style("ticks")
sns.set_context("paper")

n = 1000
k = 0.1
pc = 0.6
pd = 0.3
pn = 0.1
p = 0.1932
#G = nx.watts_strogatz_graph(n, 10, p)
G = nx.erdos_renyi_graph(n, p)
#G  = nx.barabasi_albert_graph(n, 10)
#G = nx.powerlaw_cluster_graph(1000, 10, 0.2)
A = nx.to_numpy_matrix(G)
nx.draw(G)


Nodes = np.arange(0, n, 1)
S_init = random.sample(Nodes, 995)
I_init = random.sample(np.delete(Nodes, S_init), 5)
R_init = np.delete(Nodes, np.concatenate((S_init, I_init)))

time = np.linspace(0, 100, 100)
D = {}
S = [len(S_init)]
I = [len(I_init)]
R = [len(R_init)]
a = np.zeros(20)
b = np.zeros(20)
c = np.zeros(20)
e = np.zeros(20)
Ro = np.zeros(20)
S_last_20 = np.zeros(20)
I_last_20 = np.zeros(20)
R_last_20 = np.zeros(20)

for t in range(len(time)):
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
            S_last_20[t-80] = Suscept_count
            I_last_20[t-80] = Infect_count
            R_last_20[t-80] = Recov_count
Ro_N = np.mean(Ro)

S_N = np.mean(S_last_20)
I_N = np.mean(I_last_20)
R_N = np.mean(R_last_20)
Ip = max(I)
#IIp = int(np.where(np.array(I)==Ip)[0])+1
a = np.array(I)
ind = np.unravel_index(np.argmax(a, axis = 0), a.shape)[0]

print Ro
R_N = np.mean(Ro)        
plt.figure()
plt.plot(S, color = 'b',  label = "S")
plt.plot(I, color = 'r',  label = "I")
plt.plot(R, color = 'g',  label = "R")
plt.xlabel("time", fontsize = 15)
plt.ylabel("Populations", fontsize = 15)
plt.title("SIR simulation on a SF network", fontsize = 15)
plt.tick_params(labelsize = 15)
plt.legend(prop = {"size":15})
plt.tight_layout()
plt.show()
