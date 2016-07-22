#include <iostream>
#include "ff_isom_base_change.h"
#include <flint/profiler.h>

using namespace std;

void test_base_change(slong degree) {
	cout << "degree: " << degree << "\n";
	mp_limb_t p = 9001;

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t f, g, modulus, temp;
	nmod_poly_init(f, p);
	nmod_poly_init(g, p);
	nmod_poly_init(modulus, p);
	nmod_poly_init(temp, p);

	nmod_poly_randtest(f, state, degree - 1);
	nmod_poly_randtest(temp, state, degree - 1);
	nmod_poly_randtest_monic_irreducible(modulus, state, degree);

	nmod_poly_compose_mod(g, temp, f, modulus);

	timeit_t time;
	timeit_start(time);

	FFIsomBaseChange ffIsomBaseChane;
	ffIsomBaseChane.change_basis(f, f, g, modulus);

	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	if (nmod_poly_equal(f, temp))
		cout << "ok\n";
	else
		cout << "oops\n";

	nmod_poly_clear(f);
	nmod_poly_clear(g);
	nmod_poly_clear(temp);
	nmod_poly_clear(modulus);
	flint_randclear(state);
}

int main() {

	for (slong i = 100; i < 120; i++) {
		test_base_change(i);
		cout << "----------------------\n";
	}

	return 0;
}

