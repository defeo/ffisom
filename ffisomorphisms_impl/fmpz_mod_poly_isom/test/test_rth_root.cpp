#include <iostream>
#include "cyclotomic_ext_rth_root.h"
#include "util.h"
#include <flint/profiler.h>

using namespace std;

void test_rth_root(slong degree, slong prime) {

	cout << "degree: " << degree << "\n";
	cout << "prime: " << prime << "\n";
	fmpz_t p;
	fmpz_init_set_ui(p, prime);

	flint_rand_t state;
	flint_randinit(state);

	// build the cyclotomic extension
	Util util;
	slong d = util.compute_multiplicative_order(p, degree);
	fmpz_mod_poly_t f;
	fmpz_mod_poly_init(f, p);

	for (slong i = 0; i < degree; i++)
		fmpz_mod_poly_set_coeff_ui(f, i, 1);

	fmpz_mod_poly_factor_t factors;
	fmpz_mod_poly_factor_init(factors);
	fmpz_mod_poly_factor_equal_deg(factors, f, d);

	fq_ctx_t ctx;
	fq_ctx_init_modulus(ctx, &factors->poly[0], "x");

	fq_t root;
	fq_init(root, ctx);

	fq_t rth_power;
	fq_init(rth_power, ctx);
	fq_randtest_not_zero(rth_power, state, ctx);
	fq_pow_ui(rth_power, rth_power, degree, ctx);

	timeit_t time;
	timeit_start(time);

	CyclotomicExtRthRoot cyclotomicExtRthRoot;
	cyclotomicExtRthRoot.compute_rth_root(root, rth_power, degree, ctx);

	timeit_stop(time);
	cout << "root: " << (double) time->wall / 1000.0 << "\n";

	fq_pow_ui(root, root, degree, ctx);
	if (fq_equal(root, rth_power, ctx))
		cout << "ok" << "\n";
	else
		cout << "oops" << "\n";

	fmpz_clear(p);
	flint_randclear(state);
	fmpz_mod_poly_clear(f);
	fmpz_mod_poly_factor_clear(factors);
	fq_clear(root, ctx);
	fq_clear(rth_power, ctx);
	fq_ctx_clear(ctx);
}


int main() {
	
	n_factor_t factors;
	
	for (slong i = 10; i < 30; i++) {
		slong p = n_nth_prime(i);
		n_factor_init(&factors);
		n_factor(&factors, p - 1, 1);
		test_rth_root(factors.p[factors.num - 1], p);
		cout << "\n---------------------------------------\n";
	}

	cout << "/////////////////////////////////////////////////\n";
	
	for (slong i = 10; i < 30; i++) {
		slong p = n_nth_prime(i);
		slong degree = n_nth_prime(i + 5);
		test_rth_root(degree, p);
		cout << "\n---------------------------------------\n";
	}

	return 0;
}

