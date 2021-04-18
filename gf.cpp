#include "gf.h"

GField::GField(int _poly) :
	poly(_poly)
	{
	
	m = 0;
	for (; _poly > 0; _poly >>= 1, m++);
	N = (1 << --m) - 1;

	alpha = new int[N+1];
	value = new int[N+1];
	
	for (int i = 0, x = 1; i < N+1; i++) {
		alpha[i] = x;
		value[x % N] = i;
		x <<= 1;
		if ((x & (N+1)) != 0) x ^= poly;
	}
}

GField::~GField() {
	delete[] alpha;
	delete[] value;
}

