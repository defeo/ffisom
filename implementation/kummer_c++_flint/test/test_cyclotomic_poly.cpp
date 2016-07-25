
#include <iostream>
#include "nmod_cyclotomic_poly.h"
#include "util.h"
#include <flint/profiler.h>
#include <flint/nmod_poly.h>


using namespace std;

void test_construct() {
	flint_rand_t state;
	flint_randinit(state);

	slong p = 13;
	nmod_poly_t cyclo_poly;
	nmod_poly_init(cyclo_poly, p);

	NModCyclotomicPoly nModCyclotomicPoly;

	cout << "first 30 cyclotomic polynomials with coeffs mod " << p << " :\n\n";
	for (slong i = 1; i <= 30; i++) {
		nModCyclotomicPoly.construct_cyclo(cyclo_poly, i);
		cout << i << ":  " << nmod_poly_get_str_pretty(cyclo_poly, "x") << "\n";
	}


	flint_randclear(state);
	nmod_poly_clear(cyclo_poly);
}

void test_edf() {
	flint_rand_t state;
	flint_randinit(state);

	slong p = 13;
	nmod_poly_t cyclo_poly;
	nmod_poly_init(cyclo_poly, 13);

	NModCyclotomicPoly nModCyclotomicPoly;


	for (slong i = 20; i <= 50; i += 5) {
		cout << "factors for " << i << "-th cyclotomic polynomial: \n";
		nModCyclotomicPoly.construct_cyclo(cyclo_poly, i);

		nmod_poly_factor_t factors;
		nmod_poly_factor_init(factors);
		nModCyclotomicPoly.all_irred_factors(factors, i, p);

		for (slong j = 0; j < factors->num; j++) {
			cout << j << ":  " << nmod_poly_get_str_pretty(&factors->p[j], "x") << "\n";
		}

		cout << "\n";
		nmod_poly_factor_clear(factors);
	}


	flint_randclear(state);
	nmod_poly_clear(cyclo_poly);
}

void test_single_factor() {
	flint_rand_t state;
	flint_randinit(state);

	slong p = 13;
	nmod_poly_t cyclo_poly;
	nmod_poly_init(cyclo_poly, 13);

	NModCyclotomicPoly nModCyclotomicPoly;


	for (slong i = 20; i <= 150; i += 5) {
		cout << "a factor for " << i << "-th cyclotomic polynomial: \n";
		nModCyclotomicPoly.construct_cyclo(cyclo_poly, i);
		nModCyclotomicPoly.single_irred_factor(cyclo_poly, i, p);
		cout << nmod_poly_get_str_pretty(cyclo_poly, "x") << "\n";
		cout << "\n";
	}

	flint_randclear(state);
	nmod_poly_clear(cyclo_poly);
}

int main() {

	test_single_factor();
	return 0;
}