#include <iostream>
#include "nmod_min_poly.h"
#include <flint/profiler.h>

using namespace std;

void test_minpoly(slong degree) {
	cout << "degree: " << degree << "\n";
	mp_limb_t p =  9001;

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t minpoly, f, modulus;
	nmod_poly_init(f, p);
	nmod_poly_init(minpoly, p);
	nmod_poly_init(modulus, p);
	nmod_poly_randtest(f, state, degree - 1);
	nmod_poly_randtest_monic_irreducible(modulus, state, degree);

	timeit_t time;
	timeit_start(time);
	
	NmodMinPoly nmodMinPoly;
	nmodMinPoly.minimal_polynomial(minpoly, f, modulus);
	
	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	slong irr = nmod_poly_is_irreducible(minpoly);
	nmod_poly_compose_mod(minpoly, minpoly, f, modulus);
	
	if (irr && nmod_poly_is_zero(minpoly))
		cout << "ok\n";
	else
		cout << "oops\n";

	flint_randclear(state);
	nmod_poly_clear(minpoly);
	nmod_poly_clear(f);
	nmod_poly_clear(modulus);
}


int main() {
	
	for (slong i = 100; i < 120; i++){
		test_minpoly(i);
		cout << "----------------------\n";
	}
	
	return 0;
}

