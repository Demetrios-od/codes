# Galois field and operations in it

polynomials = [0, 0x3, 0x7, 0xb, 0x13, 0x25, 0x43, 0x83, 0x11d]

def suggest_polynomial(N):
	m = 0
	while (1 << m) < N:
		m += 1
	return polynomials[m]


class gfield:

	def __init__(self, poly):
		self.poly = poly
		self.m = 0
		p = poly
		while p > 0:
			p >>= 1
			self.m += 1
		self.m -= 1
		self.N = (1 << self.m) - 1

		self.alpha = [1]*(self.N+1)
		self.value = [0]*(self.N+1)
		self.value[0] = -1

		# alpha[0] = 1
		# alpha[N] = 1
		i = 1
		x = 2
		while x != 1:
			self.alpha[i] = x
			self.value[x] = i
			i += 1
			x <<= 1
			if (x & (self.N+1) != 0):
				x ^= poly

	def powa(self, p):
		if p < 0:
			p += (-p//self.N + 1) * self.N
		return self.alpha[p % self.N]

	def log(self, x):
		return self.value[x]

	def pow(self, x, p):
		if x == 0:
			return 0
		return self.powa(self.value[x]*p)

	def mul(self, x, y):
		if x == 0 or y == 0:
			return 0
		return self.alpha[(self.value[x] + self.value[y]) % self.N]

	def inv(self, x):
		if x == 0:
			print('Inversion of zero')
			exit()
		return self.alpha[self.N - self.value[x]]

	def div(self, x, y):
		if y == 0:
			print('Division by zero')
			exit()
		if x == 0:
			return 0
		return self.alpha[(self.N + self.value[x] - self.value[y]) % self.N]

	def sqrt(self, x):
		if x == 0:
			return 0
		return self.alpha[(self.value[x] + self.N * (self.value[x] & 1)) >> 1]

	def trace(self, x):
		y = x
		z = self.value[x]
		for i in range(1, self.m):
			z = (z << 1) % self.N
			y ^= self.alpha[z]
		return y
