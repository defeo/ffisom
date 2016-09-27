#include <iostream>
#include "nmod_fast_irred.h"
#include <flint/profiler.h>

using namespace std;

void test_irred_flint(ulong p, ulong r, slong e) {

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

void test_irred_shoup(ulong p, ulong r, slong e) {
	slong n = n_pow(r, e); // overflow party
	cout << "characteristic: " << p << "\n";
	cout << "degree: " << r << "^" << e << " = " << n << "\n";

	nmod_poly_t irred;
	nmod_poly_init(irred, p);

	timeit_t time;
	timeit_start(time);

	NmodFastIrred nmodFastIrred;
	nmodFastIrred.irred_prime_power_shoup(irred, r, e, p);

	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	slong irr;

	timeit_start(time);

	irr = nmod_poly_is_irreducible_ddf(irred);

	timeit_stop(time);
	cout << "ddf: " << (double) time->wall / 1000.0 << "\n";

	/*
	timeit_start(time);

	irr = nmod_poly_is_irreducible_rabin(irred);

	timeit_stop(time);
	cout << "rabin: " << (double) time->wall / 1000.0 << "\n";
	*/

	if (irr && nmod_poly_degree(irred) == n)
		cout << "ok\n";
	else
		cout << "oops\n";

	nmod_poly_clear(irred);
}

void test_irred_adleman_lenstra(ulong p, ulong r, slong e) {
	slong n = n_pow(r, e); // overflow party
	cout << "characteristic: " << p << "\n";
	cout << "degree: " << r << "^" << e << " = " << n << "\n";

	nmod_poly_t irred;
	nmod_poly_init(irred, p);

	timeit_t time;
	timeit_start(time);

	NmodFastIrred nmodFastIrred;
	nmodFastIrred.irred_prime_power_adleman_lenstra(irred, r, e, p);

	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	slong irr;

	timeit_start(time);

	irr = nmod_poly_is_irreducible_ddf(irred);

	timeit_stop(time);
	cout << "ddf: " << (double) time->wall / 1000.0 << "\n";

	/*
	timeit_start(time);

	irr = nmod_poly_is_irreducible_rabin(irred);

	timeit_stop(time);
	cout << "rabin: " << (double) time->wall / 1000.0 << "\n";
	*/

	if (irr && nmod_poly_degree(irred) == n)
		cout << "ok\n";
	else {
		if (!irr)
			cout << "oops (irr)\n";
		else
			cout << "oops (deg)\n";
	}

	nmod_poly_clear(irred);
}

void test_irred_adleman_lenstra_factor(ulong p, ulong r, slong e) {
	slong n = n_pow(r, e); // overflow party
	cout << "characteristic: " << p << "\n";
	cout << "degree: " << r << "^" << e << " = " << n << "\n";

	nmod_poly_t irred;
	nmod_poly_init(irred, p);

	timeit_t time;
	timeit_start(time);

	NmodFastIrred nmodFastIrred;
	nmodFastIrred.irred_prime_power_adleman_lenstra_factor(irred, r, e, p);

	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	slong irr;

	timeit_start(time);

	irr = nmod_poly_is_irreducible_ddf(irred);

	timeit_stop(time);
	cout << "ddf: " << (double) time->wall / 1000.0 << "\n";

	/*
	timeit_start(time);

	irr = nmod_poly_is_irreducible_rabin(irred);

	timeit_stop(time);
	cout << "rabin: " << (double) time->wall / 1000.0 << "\n";
	*/

	if (irr && nmod_poly_degree(irred) == n)
		cout << "ok\n";
	else {
		if (!irr)
			cout << "oops (irr)\n";
		else
			cout << "oops (deg)\n";
	}

	nmod_poly_clear(irred);
}

int main() {
	for (ulong p = n_nextprime(100, 0); p < 100000; p = n_nextprime(p+1000, 0)) {
		for (ulong r = n_nextprime(12, 0); r < 50; r = n_nextprime(r, 0)) {
			if (r == p)
				continue;
			slong e = 1;
			ulong n = r;
			while(n < 5000) {
				test_irred_shoup(p, r, e);
				test_irred_adleman_lenstra_factor(p, r, e);
				//test_irred_adleman_lenstra(p, r, e);
				//test_irred_flint(p, r, e);
				cout << "----------------------\n";
				e++;
				n *= r;
			}
		}
	}

	return 0;
}

