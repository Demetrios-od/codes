# BCH encoder and decoder

import gfield


def HARD(x):
	return 1 if x < 0 else 0

def SIGN(x):
	return 1 if x>=0 else -1


def bch_generator_polynomial(gf, d):
	used_gf_elements = [0]*gf.N
	cc = []   # cyclotomic classes

	z = 0   # current class leader
	while z < d:   # for all classes, use z<gf.N, instead use z<d
		cc += [[z]]
		used_gf_elements[z] = 1
		z = (z*2) % gf.N
		while z != cc[-1][0]:
			cc[-1] += [z]
			used_gf_elements[z] = 1
			z = (z*2) % gf.N
		# find next class leader
		while True:
			z += 1
			if z >= gf.N or used_gf_elements[z] == 0:
				break
	
	# save code zeros
	poly_zeros = []
	for i in range(1, len(cc)):
		poly_zeros += cc[i]
	poly_zeros.sort()

	# make BCH generator polynomial
	poly = [1]
	for z in poly_zeros:
		poly += [0]
		a = gf.powa(z)
		for i in range(len(poly)-1, 0, -1):
			poly[i] ^= gf.mul(poly[i-1], a)
	
	# find nonzero elements in poly
	poly_powers = []
	R = len(poly)
	for i in range(R):
		if poly[i] != 0:
			poly_powers += [R-i-1]

	return {'powers': tuple(poly_powers), 'zeros': tuple(poly_zeros)}


def gauss_elim_ints(LM):
	# LM is a matrix of bits, each row is int, size [m, m+1]
	# input matrix will be changed!

	m = len(LM)
	idx = list(range(m))   # permutation indices

	for row in range(m):
		found_1 = ((LM[row] >> row)&1) == 1
		col = row
		while not found_1 and col < m:
			# put 1 into position [row, row], if possible
			row2 = row
			while row2 < m and ((LM[row2] >> col)&1) == 0:
				row2 += 1
			found_1 = row2 < m
			if found_1:
				# swap current row and found row
				bbb = LM[row]
				LM[row] = LM[row2]
				LM[row2] = bbb
			else:
				# no row with 1 in current column
				col += 1

		if not found_1:
			break

		if row < col < m:
			# swap columns
			mask = (1 << row) | (1 << col)
			for row2 in range(m):
				bbb = ((LM[row2] >> row)&1) | ((LM[row2] >> (col-1))&2)
				LM[row2] &= ~mask
				LM[row2] |= ((bbb & 1) << col) | (((bbb >> 1)&1) << row)
			# swap index
			bbb = idx[row]
			idx[row] = idx[col]
			idx[col] = bbb


		# add current row to any row up and down where 1 is at i-th position
		for row2 in range(m):
			if row2 != row and ((LM[row2] >> row)&1) == 1:
				LM[row2] ^= LM[row]
	
	return idx


def gauss_elim(LM):
	# LM is a matrix of bits, each row is list, size [m, m+1]
	# input matrix will be changed!

	m = len(LM)
	idx = list(range(m))   # permutation indices

	for row in range(m):
		found_1 = LM[row][row] == 1
		col = row
		while not found_1 and col < m:
			# put 1 into position [row, row], if possible
			row2 = row
			while row2 < m and LM[row2][col] == 0:
				row2 += 1
			found_1 = row2 < m
			if found_1:
				# swap current row and found row
				bbb = LM[row]
				LM[row] = LM[row2]
				LM[row2] = bbb
			else:
				# no row with 1 in current column
				col += 1

		if not found_1:
			break

		if row < col < m:
			# swap columns
			for row2_vct in LM:
				bbb = row2_vct[row]
				row2_vct[row] = row2_vct[col]
				row2_vct[col] = bbb
			# swap index
			bbb = idx[row]
			idx[row] = idx[col]
			idx[col] = bbb


		# add current row to any row up and down where 1 is at i-th position
		for row2 in range(m):
			if row2 != row and LM[row2][row] == 1:
				for col2 in range(row, m+1):
					LM[row2][col2] ^= LM[row][col2]
	
	return idx


class ebch:

	def __init__(self, N, T, SPC):
		self.N = N
		self.T = T
		self.SPC = SPC
		self.D = 2*self.T + 1 + self.SPC
		self.gf = gfield.gfield(gfield.suggest_polynomial(self.N - self.SPC))

		self.S = self.gf.N + self.SPC - self.N
		if (self.S < 0):
			print('BCH: code size and field size mismatch')
			exit()
		
		self.gpoly = bch_generator_polynomial(self.gf, self.D - self.SPC)
		self.R = self.gpoly['powers'][0] + self.SPC
		self.K = self.N - self.R

		# fill xi - array for t=2 decoder
		i = 0
		while self.gf.trace(self.gf.powa(i)) == 0:
			i += 1
		ak = self.gf.powa(i)

		self.xi = [0]*self.gf.m
		for i in range(self.gf.m):
			rp = self.gf.powa(i) ^ self.gf.mul(ak, self.gf.trace(self.gf.powa(i)))
			j = 0
			while (self.gf.powa(j*2) ^ self.gf.powa(j)) != rp:
				j += 1
			self.xi[i] = self.gf.powa(j)


	def encode(self, info):
		if len(info) != self.K:
			print('BCH encoder: Incorrect data size')
			exit()
		
		cw = info.copy() + [0]*self.R
		poly = self.gpoly['powers']
		F = self.R - self.SPC
		for i in range(self.K):
			if cw[i] == 1:
				for j in range(len(poly)):
					cw[i + F - poly[j]] ^= 1
		
		if self.SPC == 1:
			cw[-1] = (sum(info) + sum(cw)) & 1
		return cw[self.K:].copy()


	def calc_syndrome(self, cw, hard=False):
		hdfun = lambda x: x if hard else HARD(x)
		synd = [0] * (self.T+1)
		i = 0
		p = self.N-1 - self.SPC   # max power of received polynom
		while p >= 0:
			if hdfun(cw[i]) == 1:
				synd[0] ^= 1
				for j in range(1, self.T+1):
					synd[j] ^= self.gf.powa((2*j-1)*p)   # this works up to t=4
			i += 1
			p -= 1
		if self.SPC == 1:
			synd[0] ^= hdfun(cw[-1])
		return synd


	def update_syndrome(self, synd, pos):
		synd[0] ^= 1
		if pos < self.N - self.SPC:
			for j in range(1, self.T+1):
				synd[j] ^= self.gf.powa((2*j-1)*(self.N-1-self.SPC-pos))   # this works up to t=4
		# a bug will happen when someone try to update shortened code in a shortened position
	
	def check_syndrome(self, synd, errpos):
		synd2 = [0]*(self.T+1)
		for p in errpos:
			self.update_syndrome(synd2, p)
		return synd[1:] == synd2[1:]


	def decode_alg(self, synd):
		if all(s == 0 for s in synd[1:]):
			return True, []


		def solve4(a, b, c):
			Lx = lambda x: gpow(x,4) ^ gmul(a, gpow(x,2)) ^ gmul(b,x)

			# fill matrix LM
			m = self.gf.m
			LM = [[0]*(m+1) for i in range(m)]
			for i in range(m):
				Lai = Lx(self.gf.powa(i))
				for j in range(m):
					LM[j][i] = (Lai >> j)&1
				LM[i][m] = (c >> i)&1
			
			idx = gauss_elim(LM)

			# two last rows must be zeros, other rows must be nonzero
			if any(LM[-1]) or any(LM[-2]) or not any(LM[-3]):
				return False, []

			# find roots of the error locator polynomial
			# the system has 4 roots, but one of them was added by ourself
			roots = [0]*4
			for i in range(4):
				for j in range(m-2):
					bbb = (LM[j][m-2] & ((i>>1)&1)) ^ (LM[j][m-1] & (i&1)) ^ LM[j][m]
					roots[i] |= bbb << idx[j]
				roots[i] |= ((i>>1)&1) << idx[m-2]
				roots[i] |= (i&1) << idx[m-1]
				
			return True, roots

		
		if self.T > 4:
			print('BCH decoder: currently supported up to t=4')
			exit()
		
		gpow = lambda x,p: self.gf.pow(x,p)
		gmul = lambda x,y: self.gf.mul(x,y)
		gmul3 = lambda x,y,z: self.gf.mul(self.gf.mul(x,y), z)
		gmul4 = lambda x,y,z,w: self.gf.mul(self.gf.mul(x,y), self.gf.mul(z,w))
		gdiv = lambda x,y: self.gf.div(x,y)

		D = [0]*5
		
		if self.T > 3:
			# 4 errors
			D[0] = gmul4(synd[1], synd[2], synd[3], synd[4]) ^ \
				gmul(synd[1], gpow(synd[2], 5)) ^ \
				gmul(gpow(synd[1], 9), synd[4]) ^ \
				gmul(synd[1], gpow(synd[3], 3)) ^ \
				gmul3(gpow(synd[1], 4), synd[3], synd[4]) ^ \
				gpow(gmul(synd[1], synd[2]), 4) ^ \
				gmul3(gpow(synd[1], 8), synd[2], synd[3]) ^ \
				gpow(gmul(synd[2], synd[3]), 2) ^ \
				gmul(gpow(synd[2], 3), synd[4]) ^ \
				gpow(synd[1], 16)
			
			if D[0] != 0:
				D[1] = gmul3(gpow(synd[1], 9), synd[2], synd[3]) ^ \
					gmul3(synd[1], gpow(synd[2], 3), synd[4]) ^ \
					gmul3(gpow(synd[1], 5), synd[3], synd[4]) ^ \
					gmul(gpow(synd[1], 5), gpow(synd[2], 4)) ^ \
					gpow(synd[1], 17) ^ \
					gmul(synd[1], gpow(gmul(synd[2], synd[3]), 2)) ^ \
					gmul4(gpow(synd[1], 2), synd[2], synd[3], synd[4]) ^ \
					gmul(gpow(synd[1], 2), gpow(synd[2], 5)) ^ \
					gmul(gpow(synd[1], 10), synd[4]) ^ \
					gmul(gpow(synd[1], 2), gpow(synd[3], 3))

				D[2] = gmul(gpow(synd[1], 9), gpow(synd[2], 3)) ^ \
					gmul3(synd[1], synd[2], gpow(synd[4], 2)) ^ \
					gmul(gpow(synd[1], 13), synd[3]) ^ \
					gmul3(gpow(synd[1], 5), gpow(synd[2], 2), synd[4]) ^ \
					gmul3(synd[1], gpow(synd[2], 4), synd[3]) ^ \
					gmul3(synd[1], gpow(synd[3], 2), synd[4]) ^ \
					gmul3(gpow(synd[1], 10), synd[2], synd[3]) ^ \
					gmul3(gpow(synd[1], 2), gpow(synd[2], 3), synd[4]) ^ \
					gmul3(gpow(synd[1], 6), synd[3], synd[4]) ^ \
					gmul(gpow(synd[1], 6), gpow(synd[2], 4)) ^ \
					gpow(synd[1], 18) ^ \
					gpow(gmul3(synd[1], synd[2], synd[3]), 2) ^ \
					gmul(gpow(synd[1], 4), gpow(synd[4], 2)) ^ \
					gmul3(gpow(synd[1], 8), synd[2], synd[4]) ^ \
					gmul(synd[2], gpow(synd[3], 3)) ^ \
					gmul3(gpow(synd[2], 2), synd[3], synd[4]) ^ \
					gmul(gpow(synd[1], 8), gpow(synd[3], 2)) ^ \
					gmul(gpow(synd[1], 12), gpow(synd[2], 2))

				D[3] = gmul(gpow(synd[1], 13), gpow(synd[2], 2)) ^ \
					gmul(gpow(synd[1], 5), gpow(synd[4], 2)) ^ \
					gmul(gpow(synd[1], 9), gpow(synd[3], 2)) ^ \
					gmul(synd[1], gpow(synd[2], 6)) ^ \
					gmul(gpow(synd[1], 10), gpow(synd[2], 3)) ^ \
					gmul3(gpow(synd[1], 2), synd[2], gpow(synd[4], 2)) ^ \
					gmul(gpow(synd[1], 14), synd[3]) ^ \
					gmul3(gpow(synd[1], 6), gpow(synd[2], 2), synd[4]) ^ \
					gmul3(gpow(synd[1], 2), gpow(synd[2], 4), synd[3]) ^ \
					gmul3(gpow(synd[1], 2), gpow(synd[3], 2), synd[4]) ^ \
					gmul(gpow(synd[1], 16), synd[2]) ^ \
					gmul(gpow(synd[2], 3), gpow(synd[3], 2)) ^ \
					gmul3(gpow(synd[1], 8), gpow(synd[2], 2), synd[3]) ^ \
					gmul(gpow(synd[2], 4), synd[4]) ^ \
					gmul(gpow(synd[1], 12), synd[4]) ^ \
					gmul(gpow(synd[1], 4), gpow(synd[3], 3))

				D[4] = gmul(gpow(synd[1], 14), gpow(synd[2], 2)) ^ \
					gmul(gpow(synd[1], 6), gpow(synd[4], 2)) ^ \
					gmul(gpow(synd[1], 10), gpow(synd[3], 2)) ^ \
					gmul(gpow(synd[1], 2), gpow(synd[2], 6)) ^ \
					gpow(gmul(synd[2], synd[4]), 2) ^ \
					gpow(gmul4(synd[1], synd[1], synd[2], synd[3]), 2) ^ \
					gpow(synd[1], 20) ^ \
					gpow(synd[3], 4)

				if D[3] == 0:
					a = gdiv(D[2], D[4])
					b = gdiv(D[1], D[4])
					c = gdiv(D[0], D[4])
					res, roots = solve4(a, b, c)
				else:
					e2 = gdiv(D[1], D[3])
					e = self.gf.sqrt(e2)
					c = self.gf.inv(gmul(e2, e2 ^ gdiv(D[2], D[4])) ^ gdiv(D[0], D[4]))
					b = gmul(gdiv(D[3], D[4]), c)
					a = gmul(self.gf.sqrt(gmul(D[1], D[3])) ^ D[2], gdiv(c, D[4]))
					res, roots_y = solve4(a, b, c)
					roots = [self.gf.inv(y) ^ e for y in roots_y]
				
				pos = [self.N-1 - self.SPC - self.gf.log(self.gf.inv(rt)) for rt in roots]
				if not res or len(pos) != 4 or not self.check_syndrome(synd, pos) or any([p<0 for p in pos]):
					return False, []
				pos.sort()
				return True, pos
		
		if self.T > 2:
			# 3 errors
			D[0] = gmul3(synd[1], synd[2], synd[3]) ^ \
				gpow(synd[1], 9) ^ \
				gmul(gpow(synd[1], 4), synd[3]) ^ \
				gpow(synd[2], 3)
			
			if D[0] != 0:
				D[1] = gmul(synd[1], gpow(synd[2], 3)) ^ \
					gmul3(gpow(synd[1], 2), synd[2], synd[3]) ^ \
					gpow(synd[1], 10) ^ \
					gmul(gpow(synd[1], 5), synd[3])
				D[2] = gmul(synd[1], gpow(synd[3], 2)) ^ \
					gmul(gpow(synd[1], 8), synd[2]) ^ \
					gmul(gpow(synd[1], 2), gpow(synd[2], 3)) ^ \
					gmul(gpow(synd[2], 2), synd[3]) ^ \
					gmul(gpow(synd[1], 6), synd[3]) ^ \
					gmul(gpow(synd[1], 5), gpow(synd[2], 2))
				D[3] = gmul(gpow(synd[1], 6), gpow(synd[2], 2)) ^ \
					gmul(gpow(synd[1], 2), gpow(synd[3], 2)) ^ \
					gpow(synd[2], 4) ^ \
					gpow(synd[1], 12)
				
				L = [gdiv(x, D[0]) for x in D]

				A = self.gf.inv(L[3])
				B = gmul(A,L[1])
				C = gmul(A,L[2])   # extra root

				c = gmul(A,C)
				b = gmul(B,C) ^ A
				a = gmul(C,C) ^ B

				res, roots = solve4(a, b, c)
				pos = [self.N-1 - self.SPC - self.gf.log(self.gf.inv(rt)) for rt in roots if rt != C]
				if not res or len(pos) != 3 or not self.check_syndrome(synd, pos) or any([p<0 for p in pos]):
					return False, []
				pos.sort()
				return True, pos

		
		if self.T > 1:
			A = synd[1]
			B = gpow(synd[1], 3) ^ synd[2]
			if A != 0 and B != 0:
				u = gdiv(B, gpow(A,3))
				a = gdiv(gpow(A,2), B)
				z = 0
				for i in range(self.gf.m):
					z ^= ((u >> i)&1) * self.xi[i]
				roots = [z, z^1]
				pos = [self.N-1 - self.SPC - self.gf.log(self.gf.inv(gmul(rt,a))) for rt in roots if rt != 0]
				if len(pos) != 2 or not self.check_syndrome(synd, pos) or any([p<0 for p in pos]):
					return False, []
				pos.sort()
				return True, pos
		
		if self.T > 0 and synd[1] != 0:
			pos = [self.N-1 - self.SPC - self.gf.log(synd[1])]
			if not self.check_syndrome(synd, pos) or pos[0] < 0:
				return False, []
			return True, pos
		
		return False, []


	def decode_hard(self, synd):
		res, pos = self.decode_alg(synd)
		if res and self.SPC == 1:
			n = len(pos)
			if n < self.T and (n + synd[0]) & 1 == 1:
				pos += [self.N-1,]
			if n == self.T and (synd[0] + self.T) & 1 == 1:
				res = False
		return res, pos
