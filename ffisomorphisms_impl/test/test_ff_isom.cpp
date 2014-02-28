#include <iostream>
#include "ff_isom.h"
#include "fmpz_min_poly.h"
#include <flint/profiler.h>

using namespace std;

void test_build_isom(slong degree, slong prime) {
	cout << "extension degree: " << degree << "\n";
	cout << "characteristic: " << prime << "\n";
	fmpz_t p;
	fmpz_init_set_ui(p, prime);

	flint_rand_t state;
	flint_randinit(state);

	fmpz_mod_poly_t f1, f2;
	fmpz_mod_poly_init(f1, p);
	fmpz_mod_poly_init(f2, p);

	cout << "building finite fields...\n";
	fmpz_mod_poly_randtest_monic_irreducible(f1, state, degree + 1);
	FmpzMinPoly fmpzMinPoly;
	while (fmpz_mod_poly_degree(f2) < degree) {
		fmpz_mod_poly_randtest_not_zero(f2, state, degree);
		fmpzMinPoly.minimal_polynomial(f2, f2, f1);
	}

	cout << "building an isomorphism...\n";
	timeit_t time;
	timeit_start(time);
	FFIsomorphism ffIsomorphism(f1, f2);
	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	cout << "testing the isomorphism... ";
	fmpz_mod_poly_t temp;
	fmpz_mod_poly_init(temp, p);
	// set temp = x
	fmpz_mod_poly_set_coeff_ui(temp, 1, 1);
	// compute the image of x
	ffIsomorphism.compute_image(temp, temp);

	// the image of x should be a root of f1 in F_p[X]/<f2>
	fmpz_mod_poly_compose_mod(temp, f1, temp, f2);
	if (fmpz_mod_poly_is_zero(temp))
		cout << "ok\n";
	else
		cout << "ooops\n";

	fmpz_clear(p);
	flint_randclear(state);
	fmpz_mod_poly_clear(f1);
	fmpz_mod_poly_clear(f2);
	fmpz_mod_poly_clear(temp);
}

int main() {
	n_factor_t factors;
	for (slong i = 50; i < 60; i++) {
		slong p = n_nth_prime(i);
		n_factor_init(&factors);
		n_factor(&factors, p - 1, 1);
		test_build_isom(factors.p[factors.num - 1], p);
		cout << "\n---------------------------------------\n";
	}

	// random case
	for (slong i = 10; i < 30; i++){
		test_build_isom(n_nth_prime(i), n_nth_prime(i + 5));
		cout << "\n---------------------------------------\n";
	}

	return 0;
}

