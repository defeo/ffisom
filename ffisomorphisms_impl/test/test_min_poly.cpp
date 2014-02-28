#include <iostream>
#include "fmpz_min_poly.h"
#include <flint/profiler.h>

using namespace std;

void test_minpoly(slong degree) {
	cout << "degree: " << degree << "\n";
	fmpz_t p;
	fmpz_init_set_ui(p, 9001);

	flint_rand_t state;
	flint_randinit(state);

	fmpz_mod_poly_t minpoly, f, modulus;
	fmpz_mod_poly_init(f, p);
	fmpz_mod_poly_init(minpoly, p);
	fmpz_mod_poly_init(modulus, p);
	fmpz_mod_poly_randtest(f, state, degree - 1);
	fmpz_mod_poly_randtest_monic_irreducible(modulus, state, degree);

	timeit_t time;
	timeit_start(time);
	
	FmpzMinPoly fmpzMinPoly;
	fmpzMinPoly.minimal_polynomial(minpoly, f, modulus);
	
	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	slong irr = fmpz_mod_poly_is_irreducible(minpoly);
	fmpz_mod_poly_compose_mod(minpoly, minpoly, f, modulus);
	if (irr && fmpz_mod_poly_is_zero(minpoly))
		cout << "ok\n";
	else
		cout << "oops\n";

	fmpz_clear(p);
	flint_randclear(state);
	fmpz_mod_poly_clear(minpoly);
	fmpz_mod_poly_clear(f);
	fmpz_mod_poly_clear(modulus);
}


int main() {
	
	for (slong i = 100; i < 120; i++){
		test_minpoly(i);
		cout << "----------------------\n";
	}
	
	return 0;
}

