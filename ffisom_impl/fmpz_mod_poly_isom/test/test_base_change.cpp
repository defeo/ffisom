#include <iostream>
#include "ff_isom_base_change.h"
#include <flint/profiler.h>

using namespace std;


void test_base_change(slong degree) {
	cout << "degree: " << degree << "\n";
	fmpz_t p;
	fmpz_init_set_ui(p, 9001);

	flint_rand_t state;
	flint_randinit(state);

	fmpz_mod_poly_t f, g, modulus, temp;
	fmpz_mod_poly_init(f, p);
	fmpz_mod_poly_init(g, p);
	fmpz_mod_poly_init(modulus, p);
	fmpz_mod_poly_init(temp, p);
	
	fmpz_mod_poly_randtest(f, state, degree - 1);
	fmpz_mod_poly_randtest(temp, state, degree - 1);
	fmpz_mod_poly_randtest_monic_irreducible(modulus, state, degree);
	
	fmpz_mod_poly_compose_mod(g, temp, f, modulus);
	
	timeit_t time;
	timeit_start(time);
	
	FFIsomBaseChange ffIsomBaseChane;
	ffIsomBaseChane.change_basis(f, f, g, modulus);
	
	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";
	
	if (fmpz_mod_poly_equal(f, temp))
		cout << "ok\n";
	else
		cout << "oops\n";
	
	fmpz_clear(p);
	fmpz_mod_poly_clear(f);
	fmpz_mod_poly_clear(g);
	fmpz_mod_poly_clear(temp);
	fmpz_mod_poly_clear(modulus);
	flint_randclear(state);
}


int main() {
	
	for (slong i = 100; i < 120; i++){
		test_base_change(i);
		cout << "----------------------\n";
	}
	
	return 0;
}

