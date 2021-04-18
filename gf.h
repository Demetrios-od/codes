#pragma once

#include <exception>

class GField {
private:
	int m;   // field extension
	int N;   // number of elements - 1
	int poly;
	int *alpha;
	int *value;

public:
	GField(int _poly);
	~GField();
	
	inline int ext() { return m; }
	inline int pow(int p) {
		if (p == 0) return 1;
		if (p < 0) p += (-p/N + 1) * N;
		else if (p >= N) p %= N;
		return alpha[p];
	}
	inline int log(int a) {
		if (a <= 0) return 0;
		if (a >= N) a %= N;
		return value[a];
	}
	inline int mul(int x, int y) {
		return pow(log(x) + log(y));
	}
	inline int inv(int x) {
		return pow(-log(x));
	}
	inline int div(int x, int y) {
		return pow(log(x) - log(y));
	}
};

