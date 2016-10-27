class SymMatrix:
	def __init__(self,n):
		self.mat = []
		for i in range(n):
			self.mat.append([0 for j in range(i)])

	def set(self, i, j, val):
		if i == j:
			pass
		elif j < i:
			self.mat[i][j] = val
		else:
			self.mat[j][i] = val

	def get(self, i, j):
		if j == i:
			return 0
		elif j < i:
			return self.mat[i][j]
		else:
			return self.mat[j][i]
