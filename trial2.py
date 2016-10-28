from cvxpy import *
import numpy as np


# Esquema da topologia
# 0 1 2 3 4
# 5 6 7 8 9


# Inicialização de variáveis pelo utilizador

a_ij = float(input("Condutividade térmica entre cores? (A): "))
b = float(input("Contribuição térmica da potência do core? (B): "))
K = int(input("Horizonte temporal? (K): "))

print("Topologia da rede e temperatura inicial dos cores:")
print(" 70--70--70--95--95 ")
print(" |   |   |   |   | ")
print(" 70--70--70--95--95 ")
print("\n")
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
Tini = np.zeros(10)

for i in range(10):
	Tini[i] = 70

Tini[3] = Tini[4] = Tini[8] = Tini[9] = 95
#----------------------------------------------------------------------------

# Inicialização do vector que contém o coeficiente que dá a 'contribuição' 
# térmica da potência de cada core
B = [0.05 for i in range(10)]

#----------------------------------------------------------------------------


# Constantes do problema
PMax = 4 # potência em Watts
FMax = 1 # freq em GHz
FTarget = 0.8

#----------------------------------------------------------------------------


# Variáveis de optimização
P = Variable(10) # Vector das potências de cada core

T = [] # Matriz das temperaturas de cada core
T.append(Tini) # A temperatura inicial é fixa
for k in range(1,K): #A temperatura para cada instante k é uma variável de optimização
	T.append(Variable(10))

F = Variable(10) # Vector das frequências de cada core

#----------------------------------------------------------------------------

# Definição do problema
objective = Minimize(sum_entries(P))
constraints = []



for k in range(K-1):
	for i in range(10):
		N_temp = T[k][i] 
		for j in range(10):
			N_temp += A[i][j]*(T[k][j] - T[k][i])
		constraints.append(T[k+1][i] == (N_temp + B[i]*P[i]))

for i in range(10):
	constraints.append(P[i] >= (PMax * (square(F[i]) / square(FMax))))
	constraints.append(F[i] >= 0)
	constraints.append(F[i] <= FMax)
	constraints.append(P[i] <= PMax)


constraints.append(sum_entries(F) >= 10*FTarget)

#-----------------------------------------------------------------------------

prob = Problem(objective, constraints)
result = prob.solve()

print(prob.status)
print(prob.value)

print("T[0]: ", end="\n\t")
for i in range(5):
	print("{:.3f}".format(T[0][i]), end=" | ")
print("\n\t", end="")
for i in range(5,10):
	print("{:.3f}".format(T[0][i]), end=" |  ")

print("\n")

for i in range(1,K):
	print("T[{}]: ".format(i), end="\n\t")
	for j in range(5):
		print("{:.3f}".format(T[i][j].value), end=" | ")
	print("\n\t", end="")
	for j in range(5,10):
		print("{:.3f}".format(T[i][j].value), end=" | ")
	print("\n")

print("F: [",end=" ")
for i in range(10):
	print("{:.3f}".format(F[i].value),end=" ")
print("]")

print("P: [",end=" ")
for i in range(10):
	print("{:.3f}".format(P[i].value),end=" ")
print("]")
