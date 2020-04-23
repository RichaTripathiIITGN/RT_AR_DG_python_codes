# This code loads the data matrix and calculates the contribution index of each of
# the original features onto the eigensubspace formed by principal components. Generates a plot.
import numpy as np
from numpy import linalg as LA
from numpy.linalg import matrix_rank
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


A = np.load("Data_matrix_full_3.npy")
A = A[:, 0:6]
B = np.transpose(A)

C = np.dot(A, B)

Cov_c = np.cov(C)
rank = matrix_rank(Cov_c)
eig_vals, eig_vecs = LA.eig(Cov_c)
real_eigs = eig_vals.real
imag_eigs = eig_vals.imag
Princ_comp_matrix = np.zeros(A.shape)
for i in range(rank):
    Princ_comp_matrix[:, i] = eig_vecs[:, i]
    

Contri_index = np.zeros(6)
for i in range(6):
    S = 0
    for j in range(rank):
        S  = S + np.dot(A[:,i], Princ_comp_matrix[:,j])
    Contri_index[i] = S

##Ip = Contri_index[9]
##Iip = Contri_index[10]
##Ro = Contri_index[11]
##Contri_index[9] = Ro
##Contri_index[10] = Ip
##Contri_index[11] = Iip
x = np.array([0,1,2,3,4,5])
#fig = plt.figure()
my_xticks = [r"$cc$", r"$spl$", r"$den$", r"$dia$", r"$ave\_deg$",r"$max\_deg$"]
plt.xticks(x, my_xticks)
plt.plot(x , abs(Contri_index)/max(abs(Contri_index)), "o", markersize = 10, label = "Contribution Index")
plt.xlabel(r"Network topological and disease parameters", fontsize = 20)
plt.ylabel(r"Normalized score", fontsize = 20)
plt.title("$All \ networks$",  fontsize = 20)
plt.tick_params(labelsize = 15)
plt.yscale("log")
#plt.savefig("random_50.png")
plt.tight_layout()
#plt.legend()

'''
plt.figure()
plt.plot(A[0], "--bo", label = '0th row')
plt.plot(A[9], "--ro", label = '9th row')
plt.plot(A[19], "--ko", label = '19th row')
plt.plot(A[29], "--mo", label = '29th row')
plt.plot(A[49], "--go", label = '49th row')
plt.legend()
plt.figure()
plt.plot(real_eigs/max(real_eigs), "--bo",  label = "real part of eigen values")
plt.legend()
plt.figure()
plt.plot(imag_eigs/max(imag_eigs), "--ro", label = "imag part of eigen values")
plt.legend()
'''
plt.show()
