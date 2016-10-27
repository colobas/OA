from SymMatrix import SymMatrix

# 0 1 2 3 4
# 5 6 7 8 9

A = SymMatrix(10)

for i in range(5):
	A.set(i, i+5, 1/10)
	if i < 4:
		A.set(i, i+1, 1/10)
		A.set(i+5, i+6, 1/10)


T = [70 for i in range(10)]

T[3] = T[4] = T[8] = T[9] = 95

k = 0

def updateTemp(T, i, A):
	t_i_k = T[i]
	for j in range(10):
		t_i_k += A.get(i, j)*(T[j]-T[i])
	return t_i_k

while True:
	new_T = list(map(lambda x: updateTemp(T, x, A), range(10)))
	print("k={0} | {1}".format(k, " ".join("{:.2f}".format(x) for x in new_T)))
	T = new_T
	k += 1
