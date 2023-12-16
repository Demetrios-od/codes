#!/usr/bin/env python

import numpy as np

poly_zeros = [16,14,13,11,10,9,8,6,5,1,0]   # powers of BCH generator polynomial
N = 255
R = poly_zeros[0]
K = N-R

def generator_matrix(poly_zeros, N, K):
	R = N-K
	pos = np.array(poly_zeros)
	poly = np.zeros(R+1, int)
	poly[R-pos] = 1

	G = np.zeros([K,N], int)
	for i in range(K):
		G[i, i:i+R+1] = poly

	for i in range(K-1,1,-1):
		for j in range(i-1,-1,-1):
			if G[j,i] == 1:
				G[j,:] ^= G[i,:]
	
	return G



if __name__ == "__main__":
	G = generator_matrix(poly_zeros, N, K)
	GP = G[:,K:]
	# GP = G[list(range(1,K))+[0],K:]

	P = G.T @ np.linalg.inv((G @ G.T) % 2) @ G
	print(P.shape)



	## print generator matrix as 0-1
	#for i in range(K):
		#for j in range(R):
			#print(GP[i,j], end='')
		#print()

	# print generator matrix as ints
	for i in range(K):
		t = 0
		for j in range(R):
			if GP[i,j] == 1:
				t += 2**j
		print(t, end=', ')
	print()


	# d = np.zeros(K, int)
	# for i in range(K):
	# 	print(i, end=' ')
	# 	d[i] = 1
	# 	res = d @ GP
	# 	for j in res: print(j, end='')
	# 	print((sum(res)&1)^1)
	# 	d[i] = 0
	# print()
