#include <iostream>
#include "nmod_fast_irred.h"
#include <flint/profiler.h>

using namespace std;

void test_flint_irred(ulong p, ulong r, slong e) {

	slong n = n_pow(r, e); // overflow party

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t f;
	nmod_poly_init(f, p);

	timeit_t time;
	timeit_start(time);
	nmod_poly_randtest_monic_irreducible(f, state, n + 1);
	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	flint_randclear(state);
	nmod_poly_clear(f);
}

void test_fast_irred(ulong p, ulong r, slong e) {
	slong n = n_pow(r, e); // overflow party
	cout << "characteristic: " << p << "\n";
	cout << "degree: " << r << "^" << e << " = " << n << "\n";

	nmod_poly_t irred;
	nmod_poly_init(irred, p);

	timeit_t time;
	timeit_start(time);

	NmodFastIrred nmodFastIrred;
	nmodFastIrred.irred_prime_power(irred, r, e, p);

	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	slong irr = nmod_poly_is_irreducible(irred);

	if (irr && nmod_poly_degree(irred) == n)
		cout << "ok\n";
	else
		cout << "oops\n";

	nmod_poly_clear(irred);
}

int main() {
	for (ulong p = n_nextprime(10, 0); p < 10000; p = n_nextprime(p, 0)) {
		for (ulong r = 3; r < 1000; r = n_nextprime(r, 0)) {
			if (r == p)
				continue;
			slong e = 1;
			ulong n = r;
			while(n < 30000) {
				test_fast_irred(p, r, e);
				//test_flint_irred(p, r, e);
				cout << "----------------------\n";
				e++;
				n *= r;
			}
		}
	}

	return 0;
}

