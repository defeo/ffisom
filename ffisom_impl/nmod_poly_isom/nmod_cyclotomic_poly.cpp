
#include "nmod_cyclotomic_poly.h"
#include "util.h"
#include <flint/ulong_extras.h>

void NModCyclotomicPoly::compose(nmod_poly_t result, const nmod_poly_t f, slong n) {
	nmod_poly_t temp;
	nmod_poly_init(temp, f->mod.n);

	slong degree = nmod_poly_degree(f);
	for (slong i = 0; i <= degree; i++)
		nmod_poly_set_coeff_ui(temp, i * n, nmod_poly_get_coeff_ui(f, i));

	nmod_poly_set(result, temp);
	nmod_poly_clear(temp);
}

void NModCyclotomicPoly::construct_cyclo_prime_degree(nmod_poly_t result, slong p) {
	nmod_poly_zero(result);

	for (slong i = 0; i < p; i++)
		nmod_poly_set_coeff_ui(result, i, 1);
}

void NModCyclotomicPoly::construct_cyclo_prime_power_degree(nmod_poly_t result, slong p, slong j) {
	nmod_poly_zero(result);

	slong exponent = n_pow(p, j - 1);
	for (slong i = 0; i < p; i++)
		nmod_poly_set_coeff_ui(result, i * exponent, 1);
}

void NModCyclotomicPoly::construct_cyclo(nmod_poly_t result, slong n) {

	if (n == 1) {
		nmod_poly_zero(result);
		nmod_poly_set_coeff_ui(result, 0, result->mod.n - 1);
		nmod_poly_set_coeff_ui(result, 1, 1);

		return;
	}

	n_factor_t factors;
	n_factor_init(&factors);
	n_factor(&factors, n, 1);

	nmod_poly_t temp1;
	nmod_poly_init(temp1, result->mod.n);

	nmod_poly_t temp2;
	nmod_poly_init(temp2, result->mod.n);

	// set temp1 = x - 1
	nmod_poly_set_coeff_ui(temp1, 0, temp1->mod.n - 1);
	nmod_poly_set_coeff_ui(temp1, 1, 1);

	for (slong i = 0; i < factors.num; i++) {
		compose(temp2, temp1, factors.p[i]);
		nmod_poly_div(temp2, temp2, temp1);
		nmod_poly_set(temp1, temp2);
	}

	// this could have been in the above loop!
	// just for readability
	for (slong i = 0; i < factors.num; i++)
		n /= factors.p[i];

	compose(temp1, temp1, n);
	nmod_poly_set(result, temp1);

	nmod_poly_clear(temp1);
	nmod_poly_clear(temp2);
}

void NModCyclotomicPoly::compute_power(nmod_poly_t result, const nmod_poly_t g, slong k) {
	nmod_poly_t temp;
	nmod_poly_init(temp, g->mod.n);

	// compute p^i mod n
	slong q = n_powmod(g->mod.n, k, n);

	for (slong i = 0; i < n; i++) {
		slong b = q * i % n;
		nmod_poly_set_coeff_ui(temp, b, nmod_poly_get_coeff_ui(g, i));
	}

	nmod_poly_set(result, temp);
	nmod_poly_clear(temp);
}

void NModCyclotomicPoly::compute_trace(nmod_poly_t result, const nmod_poly_t g, slong i,
		const nmod_poly_t modulus) {

	if (i == 1) {
		nmod_poly_set(result, g);
		return;
	}

	nmod_poly_t temp1;
	nmod_poly_t temp2;
	nmod_poly_init(temp1, g->mod.n);
	nmod_poly_init(temp2, g->mod.n);

	if (i % 2 == 0) {
		compute_trace(temp1, g, i / 2, modulus);
		compute_power(temp2, temp1, i / 2);
		nmod_poly_add(temp1, temp1, temp2);
		nmod_poly_rem(result, temp1, modulus);
	} else {
		compute_trace(temp1, g, i - 1, modulus);
		compute_power(temp2, temp1, 1);
		nmod_poly_add(temp1, g, temp2);
		nmod_poly_rem(result, temp1, modulus);
	}
}

void NModCyclotomicPoly::split(nmod_poly_t f1, nmod_poly_t f2, const nmod_poly_t f,
		flint_rand_t state) {

	nmod_poly_t g1;
	nmod_poly_t g2;
	nmod_poly_init(g1, f->mod.n);
	nmod_poly_init(g2, f->mod.n);

	nmod_poly_t ONE;
	nmod_poly_init(ONE, f->mod.n);
	nmod_poly_one(ONE);

	while (true) {
		nmod_poly_randtest(g1, state, 2 * s);
		compute_trace(g1, g1, s, f);

		if (f->mod.n != 2) {
			nmod_poly_powmod_ui_binexp(g1, g1, (f->mod.n - 1) / 2, f);

			nmod_poly_gcd(g2, g1, f);
			if (nmod_poly_degree(g2) != 0 && nmod_poly_degree(g2) < nmod_poly_degree(f)) {
				nmod_poly_set(f1, g2);
				nmod_poly_div(f2, f, g2);
				break;
			}

			nmod_poly_sub(g2, g1, ONE);
			nmod_poly_gcd(g2, g2, f);
			if (nmod_poly_degree(g2) != 0 && nmod_poly_degree(g2) < nmod_poly_degree(f)) {
				nmod_poly_set(f1, g2);
				nmod_poly_div(f2, f, g2);
				break;
			}

			nmod_poly_add(g2, g1, ONE);
			nmod_poly_gcd(g2, g2, f);
			if (nmod_poly_degree(g2) != 0 && nmod_poly_degree(g2) < nmod_poly_degree(f)) {
				nmod_poly_set(f1, g2);
				nmod_poly_div(f2, f, g2);
				break;
			}
		} else {
			nmod_poly_gcd(g2, g1, f);
			if (nmod_poly_degree(g2) != 0 && nmod_poly_degree(g2) < nmod_poly_degree(f)) {
				nmod_poly_set(f1, g2);
				nmod_poly_div(f2, f, g2);
				break;
			}

			nmod_poly_add(g2, g1, ONE);
			nmod_poly_gcd(g2, g2, f);
			if (nmod_poly_degree(g2) != 0 && nmod_poly_degree(g2) < nmod_poly_degree(f)) {
				nmod_poly_set(f1, g2);
				nmod_poly_div(f2, f, g2);
				break;
			}
		}
	}

	nmod_poly_clear(g1);
	nmod_poly_clear(g2);
	nmod_poly_clear(ONE);
}

void NModCyclotomicPoly::equal_degree_fact(nmod_poly_factor_t factors, const nmod_poly_t f,
		flint_rand_t state) {
	if (nmod_poly_degree(f) == s) {
		nmod_poly_factor_insert(factors, f, 1);
		return;
	}

	nmod_poly_t f1;
	nmod_poly_t f2;
	nmod_poly_init(f1, f->mod.n);
	nmod_poly_init(f2, f->mod.n);

	split(f1, f2, f, state);
	equal_degree_fact(factors, f1, state);
	equal_degree_fact(factors, f2, state);

	nmod_poly_clear(f1);
	nmod_poly_clear(f2);
}

void NModCyclotomicPoly::single_irred_factor(nmod_poly_t factor, const nmod_poly_t f,
		flint_rand_t state) {
	if (nmod_poly_degree(f) == s) {
		nmod_poly_set(factor, f);
		return;
	}

	nmod_poly_t f1;
	nmod_poly_t f2;
	nmod_poly_init(f1, f->mod.n);
	nmod_poly_init(f2, f->mod.n);

	split(f1, f2, f, state);
	if (nmod_poly_degree(f1) < nmod_poly_degree(f2)) {
		single_irred_factor(factor, f1, state);
	} else {
		single_irred_factor(factor, f2, state);
	}

	nmod_poly_clear(f1);
	nmod_poly_clear(f2);
}

void NModCyclotomicPoly::all_irred_factors(nmod_poly_factor_t factors, slong n, slong modulus) {
	Util util;
	this->s = util.compute_multiplicative_order(modulus, n);
	this->n = n;

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t cyclo_poly;
	nmod_poly_init(cyclo_poly, modulus);

	construct_cyclo(cyclo_poly, n);
	equal_degree_fact(factors, cyclo_poly, state);

	flint_randclear(state);
	nmod_poly_clear(cyclo_poly);
}

void NModCyclotomicPoly::single_irred_factor(nmod_poly_t factor, slong n, slong modulus) {
	Util util;
	this->s = util.compute_multiplicative_order(modulus, n);
	this->n = n;

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t cyclo_poly;
	nmod_poly_init(cyclo_poly, modulus);

	construct_cyclo(cyclo_poly, n);
	single_irred_factor(factor, cyclo_poly, state);

	flint_randclear(state);
	nmod_poly_clear(cyclo_poly);
}