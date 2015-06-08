#include <iostream>
#include "cyclotomic_ext_rth_root.h"
#include "util.h"
#include <flint/profiler.h>

using namespace std;

void test_rth_root(slong degree, slong prime) {

	cout << "degree: " << degree << "\n";
	cout << "prime: " << prime << "\n";
	
	mp_limb_t p = prime;

	flint_rand_t state;
	flint_randinit(state);

	// build the cyclotomic extension
	Util util;
	slong d = util.compute_multiplicative_order(p, degree);
	nmod_poly_t f;
	nmod_poly_init(f, p);

	for (slong i = 0; i < degree; i++)
		nmod_poly_set_coeff_ui(f, i, 1);

	nmod_poly_factor_t factors;
	nmod_poly_factor_init(factors);
	nmod_poly_factor_equal_deg(factors, f, d);

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

void test_rth_root_simple(slong degree, slong prime) {

	cout << "degree: " << degree << "\n";
	cout << "prime: " << prime << "\n";
	
	mp_limb_t p = prime;

	flint_rand_t state;
	flint_randinit(state);

	mp_limb_t rth_power = n_randlimb(state) % p;
	rth_power = n_powmod(rth_power, degree, p);

	timeit_t time;
	timeit_start(time);

	CyclotomicExtRthRoot cyclotomicExtRthRoot;
	mp_limb_t root = cyclotomicExtRthRoot.compute_rth_root(rth_power, degree, p);

	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	root = n_powmod(root, degree, p);
	if (root == rth_power)
		cout << "ok" << "\n";
	else
		cout << "oops" << "\n";

	flint_randclear(state);
}

int main() {
	
//	test_rth_root_simple(13, 53);
	n_factor_t factors;
	
	for (slong i = 10; i < 30; i++) {
		slong p = n_nth_prime(i);
		n_factor_init(&factors);
		n_factor(&factors, p - 1, 1);
		test_rth_root_simple(factors.p[factors.num - 1], p);
		cout << "\n---------------------------------------\n";
	}

	cout << "/////////////////////////////////////////////////\n";
	
	for (slong i = 10; i < 30; i++) {
		slong p = n_nth_prime(i);
		slong degree = n_nth_prime(i + 5);
		test_rth_root_simple(degree, p);
		cout << "\n---------------------------------------\n";
	}

	return 0;
}

