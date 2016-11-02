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
B = np.random.uniform(low=0.01, high=0.1, size=10)

# Constantes do problema
PMax = 4 # potência em Watts
FMax = 1 # freq em GHz
FTarget = 0.8

#----------------------------------------------------------------------------


# Variáveis de optimização
P = Variable(K, 10) # Vectores das potências de cada core
T = Variable(K, 10) # Vectores das temperaturas de cada core
F = Variable(K, 10) # Vectores das freqs de cada core

#----------------------------------------------------------------------------

# Definição do problema
objective = Maximize(sum_entries(sum_entries(F, axis=1)))

constraints = []

constraints.append(T[0,range(10)] == Tini) # T inicial é fixa

for k in range(K):
        for i in range(10):
                if k > 0:
                    N_temp = T[k-1, i] 
                    for j in range(10):
                        N_temp += A[i][j]*(T[k-1, j] - T[k-1, i])

                    constraints.append(T[k, i] == (N_temp + B[i]*P[k, i]))

                constraints.append(P[k, i] >= (PMax * (square(F[k, i]) / square(FMax))))
                constraints.append(P[k, i] <= PMax)
                constraints.append(F[k, i] >= 0)
                constraints.append(F[k, i] <= FMax)

constraints.append(sum_entries(F, axis=1) >= 10*FTarget)

#-----------------------------------------------------------------------------

prob = Problem(objective, constraints)
result = prob.solve()



import random
def randomcolor():
    r = lambda: random.randint(0,255)
    return '#%02X%02X%02X' % (r(),r(),r())

fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()

Tplot = fig1.add_subplot(111)
Fplot = fig2.add_subplot(111)
Pplot = fig3.add_subplot(111)

t = range(K)
c = [randomcolor() for i in range(10)]

for i in range(10):
    s = T[range(K), i].value
    Tplot.plot(t, s, color=c[i])

for i in range(10):
    s = F[range(K), i].value
    Fplot.plot(t, s, color=c[i])

for i in range(10):
    s = P[range(K), i].value
    Pplot.plot(t, s, color=c[i])


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
