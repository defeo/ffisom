/*
 * fpz_rth_root.cpp
 *
 *  Created on: Dec 12, 2013
 *      Author: javad
 */

#include "cyclotomic_ext_rth_root.h"
#include "util.h"
#include <flint/fmpz.h>
#include <iostream>
#include <flint/profiler.h>

using namespace std;

/**
 * Computes the value $\beta_n = \alpha^p + \alpha^{p^2} + \cdots + \alpha^{p^n}$ 
 * where $\alpha \in \mathbb{F}_p[z][x]$ is such that $\beta_init = \alpha^p$.
 * After recursion level i, $\beta_i = \alpha^p + \alpha^{p^2} + \cdots + \alpha^{p^i}$
 * and $\xi_i = x^{p^i}$.
 */
void CyclotomicExtRthRoot::compute_trace(fq_nmod_poly_t beta, fq_nmod_poly_t xi, slong n, const fq_nmod_poly_t modulus) {

	if (n == 1) {
		fq_nmod_poly_set(beta, beta_init, ctx);
		fq_nmod_poly_set(xi, xi_init, ctx);
		return;
	}

	slong z_degree;
	fq_nmod_poly_t temp_beta;
	fq_nmod_poly_t temp_xi;

	fq_nmod_poly_init(temp_beta, ctx);
	fq_nmod_poly_init(temp_xi, ctx);

	if (n % 2 == 0) {

		compute_trace(beta, xi, n / 2, modulus);
		fq_nmod_poly_set(temp_beta, beta, ctx);
		fq_nmod_poly_set(temp_xi, xi, ctx);
		// compute z's degree = p^(n / 2) mod r
		z_degree = n_powmod(fmpz_fdiv_ui(fq_nmod_ctx_prime(ctx), r), n / 2, r);

	} else {

		compute_trace(beta, xi, n - 1, modulus);
		fq_nmod_poly_set(temp_beta, beta_init, ctx);
		fq_nmod_poly_set(temp_xi, xi_init, ctx);
		// compute z's degree = p mod r
		z_degree = fmpz_fdiv_ui(fq_nmod_ctx_prime(ctx), r);

	}

	compute_beta(beta, temp_xi, z_degree, modulus);
	fq_nmod_poly_add(beta, beta, temp_beta, ctx);

	compute_xi(xi, temp_xi, z_degree);

	fq_nmod_poly_clear(temp_beta, ctx);
	fq_nmod_poly_clear(temp_xi, ctx);
}

/**
 * Computes the value $trace = \alpha + \alpha^p + \alpha^{p^2} + \cdots + \alpha^{p^{s - 1}}$ 
 * where $\alpha \in \mathbb{F}_p(z)[x]$, and $s$ is the degree of the cyclotomic extension.
 */
void CyclotomicExtRthRoot::compute_trace(fq_nmod_poly_t trace, const fq_nmod_poly_t alpha, const fq_nmod_poly_t modulus) {

	if (fq_nmod_ctx_degree(ctx) == 1) {
		fq_nmod_poly_set(trace, alpha, ctx);

		return;
	}

	fq_nmod_poly_t xi;
	fq_nmod_poly_init(xi, ctx);

	compute_initials(alpha, modulus);
	compute_trace(trace, xi, fq_nmod_ctx_degree(ctx) - 1, modulus);
	fq_nmod_poly_rem(trace, trace, modulus, ctx);
	fq_nmod_poly_add(trace, trace, alpha, ctx);

	fq_nmod_poly_clear(xi, ctx);
}

/**
 * Computes $\xi_i = x^{p^i}$. It is assumed that {@code xi} and {@code old_xi}
 * are of the form $f(z)x^r'$ for some $f(z) \in \mathbb{F}_p[z]$ and $r' < r$.
 */
void CyclotomicExtRthRoot::compute_xi(fq_nmod_poly_t xi, const fq_nmod_poly_t old_xi, slong z_degree) {

	fq_nmod_t temp_coeff1;
	fq_nmod_t temp_coeff2;

	fq_nmod_init(temp_coeff1, ctx);
	fq_nmod_init(temp_coeff2, ctx);

	slong degree = fq_nmod_poly_degree(xi, ctx);
	slong old_degree = fq_nmod_poly_degree(old_xi, ctx);

	fq_nmod_poly_get_coeff(temp_coeff1, xi, degree, ctx);
	fq_nmod_poly_get_coeff(temp_coeff2, old_xi, old_degree, ctx);

	compute_simple_compose(temp_coeff1, z_degree);
	fq_nmod_pow_ui(temp_coeff2, temp_coeff2, degree, ctx);
	fq_nmod_mul(temp_coeff2, temp_coeff2, temp_coeff1, ctx);

	fq_nmod_poly_zero(xi, ctx);

	degree = degree * old_degree;

	if (degree >= r) {
		fq_nmod_pow_ui(temp_coeff1, a, degree / r, ctx);
		fq_nmod_mul(temp_coeff2, temp_coeff2, temp_coeff1, ctx);
		degree %= r;
	}

	fq_nmod_poly_set_coeff(xi, degree, temp_coeff2, ctx);

	fq_nmod_clear(temp_coeff1, ctx);
	fq_nmod_clear(temp_coeff2, ctx);
}

/**
 * Computes $\beta_i = \beta^{p^i}$. It is assumed that {@code xi} is of the form
 * $f(z)x^r'$ for some $f(z) \in \mathbb{F}_p[z]$ and $r' < r$.
 */
void CyclotomicExtRthRoot::compute_beta(fq_nmod_poly_t beta, const fq_nmod_poly_t xi, slong z_degree, const fq_nmod_poly_t modulus) {

	fq_nmod_t temp_coeff1;
	fq_nmod_t temp_coeff2;
	fq_nmod_t xi_coeff;

	fq_nmod_init(temp_coeff1, ctx);
	fq_nmod_init(temp_coeff2, ctx);
	fq_nmod_init(xi_coeff, ctx);

	slong xi_degree = fq_nmod_poly_degree(xi, ctx);
	slong beta_degree = fq_nmod_poly_degree(beta, ctx);

	compute_beta_coeffs(beta, z_degree);

	fq_nmod_poly_get_coeff(xi_coeff, xi, xi_degree, ctx);
	fq_nmod_set(temp_coeff2, xi_coeff, ctx);

	// evaluate beta at xi
	for (slong i = 1; i <= beta_degree; i++) {

		fq_nmod_poly_get_coeff(temp_coeff1, beta, i, ctx);
		fq_nmod_mul(temp_coeff1, temp_coeff1, temp_coeff2, ctx);
		fq_nmod_poly_set_coeff(beta, i, temp_coeff1, ctx);

		fq_nmod_mul(temp_coeff2, temp_coeff2, xi_coeff, ctx);
	}

	compute_simple_compose(beta, xi_degree, modulus);

	fq_nmod_clear(temp_coeff1, ctx);
	fq_nmod_clear(temp_coeff2, ctx);
	fq_nmod_clear(xi_coeff, ctx);
}

/**
 * Given $\beta = \sum_i c_i(z)x^i$, computes $\beta = \sum_i c_i(z^{z_degree})x^i$.
 * It is done by first reducing $z^{z_degree}$ to an element of the cyclotomic field,
 * and then doing modular compositions.
 */
void CyclotomicExtRthRoot::compute_beta_coeffs_small_ext(fq_nmod_poly_t beta, slong z_degree) {

	fq_nmod_t temp_coeff;
	fq_nmod_t temp_comp;

	fq_nmod_init(temp_coeff, ctx);
	fq_nmod_init(temp_comp, ctx);

	// compute z^z_degree mod g(z)
	nmod_poly_set_coeff_ui(temp_comp, z_degree, 1);
	nmod_poly_rem(temp_comp, temp_comp, ctx->modulus);

	slong beta_degree = fq_nmod_poly_degree(beta, ctx);

	for (slong i = 0; i <= beta_degree; i++) {

		fq_nmod_poly_get_coeff(temp_coeff, beta, i, ctx);
		nmod_poly_compose_mod(temp_coeff, temp_coeff, temp_comp, ctx->modulus);
		fq_nmod_poly_set_coeff(beta, i, temp_coeff, ctx);
	}

	fq_nmod_clear(temp_coeff, ctx);
	fq_nmod_clear(temp_comp, ctx);
}

/**
 * Given $\beta = \sum_i c_i(z)x^i$, computes $\beta = \sum_i c_i(z^{z_degree})x^i$.
 * It is done by computing the composition $c_i(z^{z_degree})$, and reducing it to an 
 * element of the cyclotomic field, for each coefficient.
 */
void CyclotomicExtRthRoot::compute_beta_coeffs_large_ext(fq_nmod_poly_t beta, slong z_degree) {

	fq_nmod_t temp_coeff;
	fq_nmod_init(temp_coeff, ctx);

	slong beta_degree = fq_nmod_poly_degree(beta, ctx);

	for (slong i = 0; i <= beta_degree; i++) {

		fq_nmod_poly_get_coeff(temp_coeff, beta, i, ctx);
		compute_simple_compose(temp_coeff, z_degree);
		fq_nmod_poly_set_coeff(beta, i, temp_coeff, ctx);
	}

	fq_nmod_clear(temp_coeff, ctx);
}

/**
 * Given $\beta = \sum_i c_i(z)x^i$, computes $\beta = \sum_i c_i(z^{z_degree})x^i$.
 */
void CyclotomicExtRthRoot::compute_beta_coeffs(fq_nmod_poly_t beta, slong zeta_degree) {
	Util util;
	if (util.is_small_cyclotomic_ext(r, ctx->modulus->mod.n))
		compute_beta_coeffs_small_ext(beta, zeta_degree);
	else
		compute_beta_coeffs_large_ext(beta, zeta_degree);
}

/**
 * Computes $f(z^{z_degree))$.
 */
void CyclotomicExtRthRoot::compute_simple_compose(fq_nmod_t f, slong z_degree) {

	mp_limb_t temp_coeff1 = 0;
	mp_limb_t temp_coeff2 = 0;
	fq_nmod_t result;

	fq_nmod_init(result, ctx);

	temp_coeff1 = nmod_poly_get_coeff_ui(f, 0);
	nmod_poly_set_coeff_ui(result, 0, temp_coeff1);

	slong degree = nmod_poly_degree(f);
	slong result_degree = z_degree;

	for (slong i = 1; i <= degree; i++) {

		temp_coeff1 = nmod_poly_get_coeff_ui(f, i);
		temp_coeff2 = nmod_poly_get_coeff_ui(result, result_degree);

		temp_coeff1 = nmod_add(temp_coeff1, temp_coeff2, f->mod);
		nmod_poly_set_coeff_ui(result, result_degree, temp_coeff1);

		result_degree = (result_degree + z_degree) % r;
	}

	nmod_poly_rem(result, result, ctx->modulus);
	nmod_poly_set(f, result);

	fq_nmod_clear(result, ctx);
}

/**
 * Computes $f(x^{x_degree))$.
 */
void CyclotomicExtRthRoot::compute_simple_compose(fq_nmod_poly_t f, slong x_degree, const fq_nmod_poly_t modulus) {

	fq_nmod_t temp_coeff1;
	fq_nmod_t temp_coeff2;
	fq_nmod_t a_power;
	fq_nmod_poly_t result;

	fq_nmod_init(temp_coeff1, ctx);
	fq_nmod_init(temp_coeff2, ctx);
	fq_nmod_init(a_power, ctx);
	fq_nmod_poly_init(result, ctx);

	fq_nmod_poly_get_coeff(temp_coeff1, f, 0, ctx);
	fq_nmod_poly_set_coeff(result, 0, temp_coeff1, ctx);
	fq_nmod_one(a_power, ctx);

	slong degree = fq_nmod_poly_degree(f, ctx);
	slong result_degree = x_degree;

	for (slong i = 1; i <= degree; i++) {

		fq_nmod_poly_get_coeff(temp_coeff1, f, i, ctx);

		if (result_degree >= r) {
			result_degree %= r;
			fq_nmod_mul(a_power, a_power, a, ctx);
		}

		fq_nmod_mul(temp_coeff1, temp_coeff1, a_power, ctx);
		fq_nmod_poly_get_coeff(temp_coeff2, result, result_degree, ctx);
		fq_nmod_add(temp_coeff1, temp_coeff1, temp_coeff2, ctx);
		fq_nmod_poly_set_coeff(result, result_degree, temp_coeff1, ctx);

		result_degree += x_degree;
	}

	fq_nmod_poly_set(f, result, ctx);

	fq_nmod_clear(temp_coeff1, ctx);
	fq_nmod_clear(temp_coeff2, ctx);
	fq_nmod_clear(a_power, ctx);
	fq_nmod_poly_clear(result, ctx);
}

/**
 * Compute the initial values $\beta_{init} = \alpha^p$, and $\xi_{init} = x^p$
 */
void CyclotomicExtRthRoot::compute_initials(const fq_nmod_poly_t alpha, const fq_nmod_poly_t modulus) {

	fmpz_t quotient;
	slong remainder;
	fq_nmod_t temp_coeff;

	fmpz_init(quotient);
	fq_nmod_init(temp_coeff, ctx);

	// compute xi_init = x^p
	remainder = fmpz_fdiv_ui(fq_nmod_ctx_prime(ctx), r);
	fmpz_sub_ui(quotient, fq_nmod_ctx_prime(ctx), remainder);
	fmpz_divexact_ui(quotient, quotient, r);
	fq_nmod_pow(temp_coeff, a, quotient, ctx);
	fq_nmod_poly_zero(xi_init, ctx);
	fq_nmod_poly_set_coeff(xi_init, remainder, temp_coeff, ctx);

	// compute beta_init = alpha^p
	fq_nmod_poly_set(beta_init, alpha, ctx);
	compute_beta(beta_init, xi_init, remainder, modulus);

	fmpz_clear(quotient);
	fq_nmod_clear(temp_coeff, ctx);
}

/**
 * Given a proper factor of y^r - a, computes an rth root of a. Let d = deg(factor), 
 * and let b be constant term of the factor. Then (d, r) = 1, so there are integers u, v
 * such that ud + vr = 1. This method computes (-1)^{du}b^u*a^v which is an rth root of a. 
 */
void CyclotomicExtRthRoot::compute_rth_root_from_factor(fq_nmod_t root, const fq_nmod_poly_t factor) {
	slong d = fq_nmod_poly_degree(factor, ctx);
	mp_limb_t u = 0;
	mp_limb_t v = 0;

	// compute u, v such that ud - vr = 1
	n_xgcd(&u, &v, d, r);

	fq_nmod_t temp_coeff;
	fq_nmod_t a_inv_pow;
	fq_nmod_init(temp_coeff, ctx);
	fq_nmod_init(a_inv_pow, ctx);

	fq_nmod_poly_get_coeff(temp_coeff, factor, 0, ctx);
	fq_nmod_pow_ui(temp_coeff, temp_coeff, u, ctx);
	fq_nmod_inv(a_inv_pow, a, ctx);
	fq_nmod_pow_ui(a_inv_pow, a_inv_pow, v, ctx);
	fq_nmod_mul(a_inv_pow, a_inv_pow, temp_coeff, ctx);

	if ((d % 2) && (u % 2))
		fq_nmod_neg(a_inv_pow, a_inv_pow, ctx);

	fq_nmod_set(root, a_inv_pow, ctx);

	fq_nmod_clear(temp_coeff, ctx);
	fq_nmod_clear(a_inv_pow, ctx);
}

/**
 * Given a proper factor of y^r - a, computes an rth root of a. Let d = deg(factor), 
 * and let b be the constant term of factor. Then (d, r) = 1, so there are integers u, v
 * such that ud + vr = 1. This method computes (-1)^{du}b^u*a^v which is an rth root of a. 
 */
mp_limb_t CyclotomicExtRthRoot::compute_rth_root_from_factor(const mp_limb_t a, const nmod_poly_t factor) {
	slong d = nmod_poly_degree(factor);
	mp_limb_t u = 0;
	mp_limb_t v = 0;

	// compute u, v such that ud - vr = 1
	n_xgcd(&u, &v, d, r);

	mp_limb_t b = nmod_poly_get_coeff_ui(factor, 0);
	b = nmod_pow_ui(b, u, factor->mod);
	mp_limb_t a_inv = nmod_inv(a, factor->mod);
	a_inv = nmod_pow_ui(a_inv, v, factor->mod);
	b = nmod_mul(b, a_inv, factor->mod);

	if ((d % 2) && (u % 2))
		b = factor->mod.n - b;

	return b;
}

void CyclotomicExtRthRoot::compute_rth_root(fq_nmod_t root, const fq_nmod_poly_t f) {

	fq_nmod_poly_t alpha;
	fq_nmod_poly_t trace;
	fq_nmod_poly_t f_temp;
	fmpz_t power;

	fq_nmod_poly_init(alpha, ctx);
	fq_nmod_poly_init(trace, ctx);
	fq_nmod_poly_init(f_temp, ctx);
	fq_nmod_poly_set(f_temp, f, ctx);

	flint_rand_t state;
	flint_randinit(state);

	while (true) {

		fq_nmod_poly_zero(trace, ctx);

		// find a nonzero trace
		while (fq_nmod_poly_is_zero(trace, ctx)) {
			fq_nmod_poly_randtest(alpha, state, fq_nmod_poly_degree(f_temp, ctx), ctx);
			compute_trace(trace, alpha, f_temp);
		}

		// compute power = (p - 1) / 2
		fmpz_sub_ui(power, fq_nmod_ctx_prime(ctx), 1);
		fmpz_divexact_ui(power, power, 2);

		// compute trace = trace^power
		fq_nmod_poly_powmod_fmpz_binexp(trace, trace, power, f_temp, ctx);

		fq_nmod_poly_one(alpha, ctx);
		fq_nmod_poly_sub(trace, trace, alpha, ctx);
		fq_nmod_poly_gcd(alpha, trace, f_temp, ctx);

		slong degree1 = fq_nmod_poly_degree(alpha, ctx);
		slong degree2 = fq_nmod_poly_degree(f_temp, ctx);
		// check if alpha is a proper factor of f_temp
		if (0 < degree1 && degree1 < degree2) {
			if (n_gcd_full(degree1, r) != 1) {
				fq_nmod_poly_set(f_temp, alpha, ctx);
			} else {
				compute_rth_root_from_factor(root, alpha);
				break;
			}
		}
	}

	fq_nmod_poly_clear(alpha, ctx);
	fq_nmod_poly_clear(trace, ctx);
	fq_nmod_poly_clear(f_temp, ctx);
	fmpz_clear(power);
	flint_randclear(state);
}

mp_limb_t CyclotomicExtRthRoot::compute_rth_root(const nmod_poly_t f) {

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t a;
	nmod_poly_t b;
	nmod_poly_t f_temp;

	mp_limb_t root;

	nmod_poly_init(a, f->mod.n);
	nmod_poly_init(b, f->mod.n);
	nmod_poly_init(f_temp, f->mod.n);
	nmod_poly_set(f_temp, f);

	while (true) {

		nmod_poly_zero(a);

		// choose a random a that is not in F_p
		while (nmod_poly_degree(a) < 1)
			nmod_poly_randtest(a, state, nmod_poly_degree(f_temp));

		nmod_poly_gcd(b, a, f_temp);
		slong degree = nmod_poly_degree(b);
		if (!nmod_poly_is_one(b)) {
			if (n_gcd_full(degree, r) != 1) {
				nmod_poly_set(f_temp, b);
			} else {
				root = compute_rth_root_from_factor(f_temp->mod.n - nmod_poly_get_coeff_ui(f, 0), b);
				break;
			}
		}

		nmod_poly_powmod_ui_binexp(b, a, (f_temp->mod.n - 1) / 2, f_temp);

		// compute b - 1
		mp_limb_t c = nmod_poly_get_coeff_ui(b, 0);
		if (c == 0)
			c = f->mod.n - 1;
		else
			c--;
		nmod_poly_set_coeff_ui(b, 0, c);

		nmod_poly_gcd(b, b, f_temp);
		degree = nmod_poly_degree(b);
		if (!nmod_poly_is_one(b) && !nmod_poly_equal(b, f_temp)) {
			if (n_gcd_full(degree, r) != 1) {
				nmod_poly_set(f_temp, b);
			} else {
				root = compute_rth_root_from_factor(f_temp->mod.n - nmod_poly_get_coeff_ui(f, 0), b);
				break;
			}
		}
	}

	flint_randclear(state);
	nmod_poly_clear(a);
	nmod_poly_clear(b);

	return root;
}

/**
 * Supports aliasing.
 */
void CyclotomicExtRthRoot::compute_rth_root(fq_nmod_t root, const fq_nmod_t a, slong r, const fq_nmod_ctx_t ctx) {

	this->ctx = ctx;
	this->r = r;
	fq_nmod_init(this->a, ctx);
	fq_nmod_set(this->a, a, ctx);
	fq_nmod_poly_init(beta_init, ctx);
	fq_nmod_poly_init(xi_init, ctx);

	fq_nmod_poly_t f;
	fq_nmod_t temp_coeff;

	fq_nmod_poly_init(f, ctx);
	fq_nmod_init(temp_coeff, ctx);

	fq_nmod_one(temp_coeff, ctx);
	fq_nmod_poly_set_coeff(f, r, temp_coeff, ctx);
	fq_nmod_neg(temp_coeff, a, ctx);
	fq_nmod_poly_set_coeff(f, 0, temp_coeff, ctx);

	compute_rth_root(root, f);

	fq_nmod_poly_clear(f, ctx);
	fq_nmod_clear(temp_coeff, ctx);
	fq_nmod_clear(this->a, ctx);
	fq_nmod_poly_clear(beta_init, ctx);
	fq_nmod_poly_clear(xi_init, ctx);
	this->ctx = NULL;
}

mp_limb_t CyclotomicExtRthRoot::compute_rth_root(const mp_limb_t c, slong r, slong p) {
	this->r = r;
	nmod_poly_t f;
	nmod_poly_init(f, p);

	nmod_poly_set_coeff_ui(f, 0, p - c);
	nmod_poly_set_coeff_ui(f, r, 1);

	mp_limb_t root = compute_rth_root(f);

	nmod_poly_clear(f);

	return root;
}
