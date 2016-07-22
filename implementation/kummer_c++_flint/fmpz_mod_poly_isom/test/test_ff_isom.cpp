#include <iostream>
#include "ff_isom.h"
#include "util.h"
#include "fmpz_min_poly.h"
#include <flint/profiler.h>

using namespace std;

void test_build_isom(slong degree, slong prime) {

	Util util;
	cout << "characteristic: " << prime << "\n";
	cout << "extension degree: " << degree << "\n";
	cout << "cyclotomic extension degree: " << util.compute_multiplicative_order(prime, degree) << "\n";

	fmpz_t p;
	fmpz_init_set_ui(p, prime);

	flint_rand_t state;
	flint_randinit(state);

	fmpz_mod_poly_t f1, f2;
	fmpz_mod_poly_init(f1, p);
	fmpz_mod_poly_init(f2, p);

	cout << "building finite fields...\n";
	fmpz_mod_poly_randtest_monic_irreducible(f1, state, degree + 1);
	NmodMinPoly fmpzMinPoly;
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

double test_build_isom1(slong degree, slong prime) {

	fmpz_t p;
	fmpz_init_set_ui(p, prime);

	flint_rand_t state;
	flint_randinit(state);

	fmpz_mod_poly_t f1, f2;
	fmpz_mod_poly_init(f1, p);
	fmpz_mod_poly_init(f2, p);

	fmpz_mod_poly_randtest_monic_irreducible(f1, state, degree + 1);
	NmodMinPoly fmpzMinPoly;
	while (fmpz_mod_poly_degree(f2) < degree) {
		fmpz_mod_poly_randtest_not_zero(f2, state, degree);
		fmpzMinPoly.minimal_polynomial(f2, f2, f1);
	}

	timeit_t time;
	timeit_start(time);
	FFIsomorphism ffIsomorphism(f1, f2);
	timeit_stop(time);

	fmpz_clear(p);
	flint_randclear(state);
	fmpz_mod_poly_clear(f1);
	fmpz_mod_poly_clear(f2);

	return (double) time->wall / 1000.0;
}

int main() {

	flint_rand_t state;
	flint_randinit(state);

	cout << "\n/////////////////////// random case //////////////////////////\n\n";
	cout << "prime  extDegree  time \n";	
	for (slong i = 30; i < 60; i++) {
		slong prime = n_nth_prime(5 + n_randint(state, i));
		slong degree = n_nth_prime(5 + n_randint(state, i));
		double time = test_build_isom1(degree, prime);
		cout << prime << "  " << degree << "  " << time << "\n";
	}
	flint_randclear(state);
	
	cout << "\n/////////////////////// no cyclotomic extension //////////////////////////\n\n";
	cout << "prime  extDegree  time \n";
	n_factor_t factors;
	for (slong i = 100; i < 120; i++) {
		slong prime = n_nth_prime(i);
		n_factor_init(&factors);
		n_factor(&factors, prime - 1, 1);
		slong degree = factors.p[factors.num - 1];
		double time = test_build_isom1(degree, prime);
		cout << prime << "  " << degree << "  " << time << "\n";
	}

	Util util;
	cout << "\n/////////////////////// large cyclotomic extension //////////////////////////\n\n";
	cout << "prime  extDegree  time \n";
	for (slong i = 20; i < 50; i++) {
		slong prime = n_nth_prime(i);
		slong degree = n_nth_prime(i - 5);
		if (util.compute_multiplicative_order(prime, degree) == degree - 1) {
			double time = test_build_isom1(degree, prime);
			cout << prime << "  " << degree << "  " << time << "\n";
		}
	}

	return 0;
}

