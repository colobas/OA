from cvxpy import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import settings2 as st


K=250
tot_cores=10
N=5
M=2

a_ij = 1/tot_cores

# Inicialização da matriz que contém os coeficientes de condutividade térmica.
# Estes são zero para os elementos que não são adjacentes.

A = np.zeros((tot_cores, tot_cores), dtype=np.float)

for i in range(N):
    for j in range(M):
        idx = i*M + j
        if idx < (N-1)*M:
            A[idx][idx+M] = A[idx+M][idx] = a_ij
        if j < M-1:
            A[idx][idx+1] = A[idx+1][idx] = a_ij


# ---------------------------------------------------------------------------
# Constantes do problema
TMax = 100
PMax = 4 # potência em Watts
FMax = 1 # freq em GHz


# Inicialização do vector de temperaturas iniciais dos cores
Tini = st.Tini
B1 = st.B1
B2 = st.B2
B3 = st.B3

#----------------------------------------------------------------------------


# Variáveis de optimização
P1 = Variable(K, tot_cores) # Vectores das potências de cada core
T1 = Variable(K, tot_cores) # Vectores das temperaturas de cada core
P2 = Variable(K, tot_cores) # Vectores das potências de cada core
T2 = Variable(K, tot_cores) # Vectores das temperaturas de cada core
P3 = Variable(K, tot_cores) # Vectores das potências de cada core
T3 = Variable(K, tot_cores) # Vectores das temperaturas de cada core
F = Variable(K, tot_cores) # Vectores das freqs de cada core
#---------------------------------------------------------------


print("A correr trial sem influência de uma fonte externa de calor para 3 conjuntos de cores")

# Definição do problema
objective = Maximize(sum_entries(F))

constraints = []

constraints.append(T1[0,range(tot_cores)] == Tini) # T inicial é fixa
constraints.append(T2[0,range(tot_cores)] == Tini) # T inicial é fixa
constraints.append(T3[0,range(tot_cores)] == Tini) # T inicial é fixa

for k in range(K):
    if k > 0:
        for i in range(tot_cores):
            N_temp1 = T1[k-1, i]
            N_temp2 = T2[k-1, i]
            N_temp3 = T3[k-1, i]
            for j in range(tot_cores):
                N_temp1 += A[i][j]*(T1[k-1, j] - T1[k-1, i])
                N_temp2 += A[i][j]*(T2[k-1, j] - T2[k-1, i])
                N_temp3 += A[i][j]*(T3[k-1, j] - T3[k-1, i])


            constraints.append(T1[k, i] == (N_temp1 + B1[i]*P1[k, i]))
            constraints.append(T2[k, i] == (N_temp2 + B2[i]*P2[k, i]))
            constraints.append(T3[k, i] == (N_temp3 + B3[i]*P3[k, i]))
    # -------------
    constraints.append(P1[k,:] >= (PMax * (square(F[k,:] / FMax))))
    constraints.append(P2[k,:] >= (PMax * (square(F[k,:] / FMax))))
    constraints.append(P3[k,:] >= (PMax * (square(F[k,:] / FMax))))
    constraints.append(P1[k,:] <= PMax)
    constraints.append(P2[k,:] <= PMax)
    constraints.append(P3[k,:] <= PMax)
    constraints.append(F[k,:] >= 0)
    constraints.append(F[k,:] <= FMax)
    constraints.append(T1[k,:] <= TMax)
    constraints.append(T2[k,:] <= TMax)
    constraints.append(T3[k,:] <= TMax)


#-----------------------------------------------------------------------------

prob = Problem(objective, constraints)
result = prob.solve()
print(result)

import random
def randomcolor():
    r = lambda: random.randint(0,255)
    return '#%02X%02X%02X' % (r(),r(),r())




plt.figure(1)

c = [randomcolor() for i in range(tot_cores)]


plt.subplot(711)
plt.title('T1(k)')
for i in range(tot_cores):
    s = T1[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(712)
plt.title('T2(k)')
for i in range(tot_cores):
    s = T2[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(713)
plt.title('T3(k)')
for i in range(tot_cores):
    s = T3[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(714)
plt.title('F(k)')
for i in range(tot_cores):
    s = F[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(715)
plt.title('P1(k)')
for i in range(tot_cores):
    s = P1[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(716)
plt.title('P2(k)')
for i in range(tot_cores):
    s = P2[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(717)
plt.title('P3(k)')
for i in range(tot_cores):
    s = P3[range(K), i].value
    plt.plot(s, color=c[i])

plt.show()
