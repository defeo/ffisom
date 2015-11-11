#include <iostream>
#include "cyclotomic_ext_rth_root.h"
#include "util.h"
#include "nmod_cyclotomic_poly.h"
#include <flint/profiler.h>

using namespace std;

void test_rth_root(slong v, slong d, slong prime) {

	cout << "degree: " << v << "^" << d << "\n";
	cout << "prime: " << prime << "\n";

	slong degree = n_pow(v, d);
	mp_limb_t p = prime;

	flint_rand_t state;
	flint_randinit(state);

	// build the cyclotomic extension
	Util util;
	ulong s = util.compute_multiplicative_order(p, degree);
	nmod_poly_t f;
	nmod_poly_init(f, p);

	NModCyclotomicPoly nModCyclotomicPoly;
	nModCyclotomicPoly.construct_cyclo_prime_power_degree(f, v, d);
	nmod_poly_factor_t factors;
	nmod_poly_factor_init(factors);
	nmod_poly_factor_equal_deg(factors, f, s);

	fq_nmod_ctx_t ctx;
	fq_nmod_ctx_init_modulus(ctx, &factors->p[0], "x");

	fq_nmod_t root;
	fq_nmod_init(root, ctx);

	fq_nmod_t rth_power;
	fq_nmod_init(rth_power, ctx);
	fq_nmod_randtest_not_zero(rth_power, state, ctx);
	fq_nmod_pow_ui(rth_power, rth_power, degree, ctx);

	timeit_t time;
	timeit_start(time);

	CyclotomicExtRthRoot cyclotomicExtRthRoot;
	cyclotomicExtRthRoot.compute_rth_root(root, rth_power, degree, ctx);

	timeit_stop(time);
	cout << "root: " << (double) time->wall / 1000.0 << "\n";

	fq_nmod_pow_ui(root, root, degree, ctx);
	if (fq_nmod_equal(root, rth_power, ctx))
		cout << "ok" << "\n";
	else
		cout << "oops" << "\n";

	flint_randclear(state);
	nmod_poly_clear(f);
	nmod_poly_factor_clear(factors);
	fq_nmod_clear(root, ctx);
	fq_nmod_clear(rth_power, ctx);
	fq_nmod_ctx_clear(ctx);
}

int main() {

	cout << "/////////////////////////////////////////////////\n";

	for (slong i = 10; i < 30; i++) {
		slong p = n_nth_prime(i);
		slong v = n_nth_prime(i + 5);
		slong d = 1;
		test_rth_root(v, d, p);
		cout << "\n---------------------------------------\n";
	}

	cout << "/////////////////////////////////////////////////\n";

	for (slong i = 2; i < 8; i++) {
		slong p = n_nth_prime(i);
		slong v = n_nth_prime(i + 1);
		slong d = 2;
		test_rth_root(v, d, p);
		cout << "\n---------------------------------------\n";
	}

	cout << "/////////////////////////////////////////////////\n";

	for (slong i = 2; i < 5; i++) {
		slong p = n_nth_prime(i);
		slong v = n_nth_prime(i + 1);
		slong d = 3;
		test_rth_root(v, d, p);
		cout << "\n---------------------------------------\n";
	}

	return 0;
}

