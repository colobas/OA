from cvxpy import *
import numpy as np
import matplotlib.pyplot as plt

# Esquema da topologia
# 0 1 2 3 4
# 5 6 7 8 9


# Inicialização de variáveis pelo utilizador

a_ij = float(input("Condutividade térmica entre cores? (A): "))
K = int(input("Horizonte temporal? (K): "))

a = input("Prima qualquer tecla para continuar")

# Inicialização da matriz que contém os coeficientes de condutividade térmica.
# Estes são zero para os elementos que não são adjacentes.

A = np.zeros((10,10), dtype=np.float)

for i in range(5):
        A[i][i+5] = A[i+5][i] = a_ij
        if i < 4:
                A[i][i+1] = A[i+1][i] = a_ij
                A[i+5][i+6] = A[i+6][i+5] = a_ij

# ---------------------------------------------------------------------------

# Inicialização do vector de temperaturas iniciais dos cores
Tini = np.random.uniform(low=70.0, high=100.0, size=10)
B = np.random.uniform(low=0.07, high=0.1, size=10)

# Constantes do problema
TMax = 100
PMax = 4 # potência em Watts
FMax = 1 # freq em GHz

#----------------------------------------------------------------------------


# Variáveis de optimização
P = Variable(K, 10) # Vectores das potências de cada core
T = Variable(K, 10) # Vectores das temperaturas de cada core
F = Variable(K, 10) # Vectores das freqs de cada core
FTarget = Variable(1)
#----------------------------------------------------------------------------

# Definição do problema
objective = Maximize(FTarget)

constraints = []

constraints.append(T[0,range(10)] == Tini) # T inicial é fixa

for k in range(K):
    if k > 0:
        for i in range(10):
            N_temp = T[k-1, i] 
            for j in range(10):
                N_temp += A[i][j]*(T[k-1, j] - T[k-1, i])

            constraints.append(T[k, i] == (N_temp + B[i]*P[k, i]))

    constraints.append(P[k,:] >= (PMax * (square(F[k,:] / FMax))))
    constraints.append(P[k,:] <= PMax)
    constraints.append(F[k,:] >= 0)
    constraints.append(F[k,:] <= FMax)
    constraints.append(T[k,:] <= TMax)
    constraints.append(sum_entries(F[k, :]) >= 10*FTarget)


#-----------------------------------------------------------------------------

prob = Problem(objective, constraints)
result = prob.solve()
print(result)

import random
def randomcolor():
    r = lambda: random.randint(0,255)
    return '#%02X%02X%02X' % (r(),r(),r())




plt.figure(1)

c = [randomcolor() for i in range(10)]


plt.subplot(311)
plt.title('T(k)')
for i in range(10):
    s = T[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(312)
plt.title('F(k)')
for i in range(10):
    s = F[range(K), i].value
    plt.plot(s, color=c[i])

plt.subplot(313)
plt.title('P(k)')
for i in range(10):
    s = P[range(K), i].value
    plt.plot(s, color=c[i])


plt.show()






# Pretty printing daqui para a frente


#for i in range(K):
#        print("K:{} ".format(i), end="\n\t")
#        for j in range(5):
#            print("T:{:.3f} F:{:.3f} P:{:.3f}".format(T[i, j].value, F[i, j].value, P[i, j].value), end=" | ")
#        print("\n\t", end="")
#        for j in range(5,10):
#            print("T:{:.3f} F:{:.3f} P:{:.3f}".format(T[i, j].value, F[i, j].value, P[i, j].value), end=" | ")
#        print("\n")
