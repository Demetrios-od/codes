#include <cstdio>

#include "gf.h"
//#include "bch.h"

using namespace std;

int main(int argc, char **argv) {
	GField gf(0x13);
	int m = gf.ext();
	int N = 1 << m;
	printf("  n  a^n log inv\n");
	for (int i = -3; i < N*2+2; i++) {
		printf("%3d| %3d %3d %3d %3d\n", i, gf.pow(i), gf.log(i), gf.inv(i), gf.mul(i, gf.inv(i)));
	}
	
	return 0;
}
