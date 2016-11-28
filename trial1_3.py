from cvxpy import *
import numpy as np
import matplotlib.pyplot as plt
import sys
import settings1 as st


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
B = st.B
p = 0.1
print("Correndo trial para o qual há influência de uma fonte de calor externa e se tenta minimizar"
      " a diferenca de temperatura entre cada dois instantes de tempo")




Tmb = np.zeros(K)
for i in range(round(K/3), round(2*K/3)):
    Tmb[i] = 10


# Variáveis de optimização
P = Variable(K, tot_cores) # Vectores das potências de cada core
T = Variable(K, tot_cores) # Vectores das temperaturas de cada core
F = Variable(K, tot_cores) # Vectores das freqs de cada core
SumF = Variable(1)
OscTerm = Variable(1)
#----------------------------------------------------------------------------

# Definição do problema
objective = Minimize(-SumF + tv(F.T)*p)

constraints = []

constraints.append(T[0,range(tot_cores)] == Tini) # T inicial é fixa
constraints.append(SumF == sum_entries(F))
#constraints.append(OscTerm == tv(F.T))

for k in range(K):
    if k > 0:
        for i in range(tot_cores):
            N_temp = T[k-1, i] 
            for j in range(tot_cores):
                N_temp += A[i][j]*(T[k-1, j] - T[k-1, i])

            constraints.append(T[k, i] == (N_temp + B[i]*P[k, i]) + a_ij*Tmb[k-1]/10)

    constraints.append(P[k,:] >= (PMax * (square(F[k,:] / FMax))))
    constraints.append(P[k,:] <= PMax)
    constraints.append(F[k,:] >= 0)
    constraints.append(F[k,:] <= FMax)
    constraints.append(T[k,:] <= TMax)


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


plt.subplot(311)
plt.title('T(k)')
for i in range(tot_cores):
    s = T[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(312)
plt.title('F(k)')
for i in range(tot_cores):
    s = F[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(313)
plt.title('P(k)')
for i in range(tot_cores):
    s = P[range(K), i].value
    plt.plot(s, color=c[i])


plt.show()
