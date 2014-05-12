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
void CyclotomicExtRthRoot::compute_trace(fq_poly_t beta, fq_poly_t xi, slong n, const fq_poly_t modulus) {

	if (n == 1) {
		fq_poly_set(beta, beta_init, ctx);
		fq_poly_set(xi, xi_init, ctx);
		return;
	}

	slong z_degree;
	fq_poly_t temp_beta;
	fq_poly_t temp_xi;

	fq_poly_init(temp_beta, ctx);
	fq_poly_init(temp_xi, ctx);

	if (n % 2 == 0) {

		compute_trace(beta, xi, n / 2, modulus);
		fq_poly_set(temp_beta, beta, ctx);
		fq_poly_set(temp_xi, xi, ctx);
		// compute z's degree = p^(n / 2) mod r
		z_degree = n_powmod(fmpz_fdiv_ui(fq_ctx_prime(ctx), r), n / 2, r);

	} else {

		compute_trace(beta, xi, n - 1, modulus);
		fq_poly_set(temp_beta, beta_init, ctx);
		fq_poly_set(temp_xi, xi_init, ctx);
		// compute z's degree = p mod r
		z_degree = fmpz_fdiv_ui(fq_ctx_prime(ctx), r);

	}

	compute_beta(beta, temp_xi, z_degree, modulus);
	fq_poly_add(beta, beta, temp_beta, ctx);

	compute_xi(xi, temp_xi, z_degree);

	fq_poly_clear(temp_beta, ctx);
	fq_poly_clear(temp_xi, ctx);
}

/**
 * Computes the value $trace = \alpha + \alpha^p + \alpha^{p^2} + \cdots + \alpha^{p^{s - 1}}$ 
 * where $\alpha \in \mathbb{F}_p(z)[x]$, and $s$ is the degree of the cyclotomic extension.
 */
void CyclotomicExtRthRoot::compute_trace(fq_poly_t trace, const fq_poly_t alpha, const fq_poly_t modulus) {

	if (fq_ctx_degree(ctx) == 1) {
		fq_poly_set(trace, alpha, ctx);

		return;
	}

	fq_poly_t xi;
	fq_poly_init(xi, ctx);

	compute_initials(alpha, modulus);
	compute_trace(trace, xi, fq_ctx_degree(ctx) - 1, modulus);
	fq_poly_rem(trace, trace, modulus, ctx);
	fq_poly_add(trace, trace, alpha, ctx);

	fq_poly_clear(xi, ctx);
}

/**
 * Computes $\xi_i = x^{p^i}$. It is assumed that {@code xi} and {@code old_xi}
 * are of the form $f(z)x^r'$ for some $f(z) \in \mathbb{F}_p[z]$ and $r' < r$.
 */
void CyclotomicExtRthRoot::compute_xi(fq_poly_t xi, const fq_poly_t old_xi, slong z_degree) {

	fq_t temp_coeff1;
	fq_t temp_coeff2;

	fq_init(temp_coeff1, ctx);
	fq_init(temp_coeff2, ctx);

	slong degree = fq_poly_degree(xi, ctx);
	slong old_degree = fq_poly_degree(old_xi, ctx);

	fq_poly_get_coeff(temp_coeff1, xi, degree, ctx);
	fq_poly_get_coeff(temp_coeff2, old_xi, old_degree, ctx);

	compute_simple_compose(temp_coeff1, z_degree);
	fq_pow_ui(temp_coeff2, temp_coeff2, degree, ctx);
	fq_mul(temp_coeff2, temp_coeff2, temp_coeff1, ctx);

	fq_poly_zero(xi, ctx);

	degree = degree * old_degree;

	if (degree >= r) {
		fq_pow_ui(temp_coeff1, a, degree / r, ctx);
		fq_mul(temp_coeff2, temp_coeff2, temp_coeff1, ctx);
		degree %= r;
	}

	fq_poly_set_coeff(xi, degree, temp_coeff2, ctx);

	fq_clear(temp_coeff1, ctx);
	fq_clear(temp_coeff2, ctx);
}

/**
 * Computes $\beta_i = \beta^{p^i}$. It is assumed that {@code xi} is of the form
 * $f(z)x^r'$ for some $f(z) \in \mathbb{F}_p[z]$ and $r' < r$.
 */
void CyclotomicExtRthRoot::compute_beta(fq_poly_t beta, const fq_poly_t xi, slong z_degree, const fq_poly_t modulus) {

	fq_t temp_coeff1;
	fq_t temp_coeff2;
	fq_t xi_coeff;

	fq_init(temp_coeff1, ctx);
	fq_init(temp_coeff2, ctx);
	fq_init(xi_coeff, ctx);

	slong xi_degree = fq_poly_degree(xi, ctx);
	slong beta_degree = fq_poly_degree(beta, ctx);

	compute_beta_coeffs(beta, z_degree);

	fq_poly_get_coeff(xi_coeff, xi, xi_degree, ctx);
	fq_set(temp_coeff2, xi_coeff, ctx);

	// evaluate beta at xi
	for (slong i = 1; i <= beta_degree; i++) {

		fq_poly_get_coeff(temp_coeff1, beta, i, ctx);
		fq_mul(temp_coeff1, temp_coeff1, temp_coeff2, ctx);
		fq_poly_set_coeff(beta, i, temp_coeff1, ctx);

		fq_mul(temp_coeff2, temp_coeff2, xi_coeff, ctx);
	}

	compute_simple_compose(beta, xi_degree, modulus);

	fq_clear(temp_coeff1, ctx);
	fq_clear(temp_coeff2, ctx);
	fq_clear(xi_coeff, ctx);
}

/**
 * Given $\beta = \sum_i c_i(z)x^i$, computes $\beta = \sum_i c_i(z^{z_degree})x^i$.
 * It is done by first reducing $z^{z_degree}$ to an element of the cyclotomic field,
 * and then doing modular compositions.
 */
void CyclotomicExtRthRoot::compute_beta_coeffs_small_ext(fq_poly_t beta, slong z_degree) {

	fq_t temp_coeff;
	fmpz_mod_poly_t temp_comp1;
	fmpz_mod_poly_t temp_comp2;

	fq_init(temp_coeff, ctx);
	fmpz_mod_poly_init(temp_comp1, fq_ctx_prime(ctx));
	fmpz_mod_poly_init(temp_comp2, fq_ctx_prime(ctx));

	// compute z^z_degree mod g(z)
	fmpz_mod_poly_set_coeff_ui(temp_comp1, z_degree, 1);
	fmpz_mod_poly_rem(temp_comp1, temp_comp1, ctx->modulus);

	slong beta_degree = fq_poly_degree(beta, ctx);

	for (slong i = 0; i <= beta_degree; i++) {

		fq_poly_get_coeff(temp_coeff, beta, i, ctx);
		fmpz_mod_poly_set_fmpz_poly(temp_comp2, temp_coeff);
		fmpz_mod_poly_compose_mod(temp_comp2, temp_comp2, temp_comp1, ctx->modulus);
		fmpz_mod_poly_get_fmpz_poly(temp_coeff, temp_comp2);
		fq_poly_set_coeff(beta, i, temp_coeff, ctx);
	}

	fq_clear(temp_coeff, ctx);
	fmpz_mod_poly_clear(temp_comp1);
	fmpz_mod_poly_clear(temp_comp2);
}

/**
 * Given $\beta = \sum_i c_i(z)x^i$, computes $\beta = \sum_i c_i(z^{z_degree})x^i$.
 * It is done by computing the composition $c_i(z^{z_degree})$, and reducing it to an 
 * element of the cyclotomic field, for each coefficient.
 */
void CyclotomicExtRthRoot::compute_beta_coeffs_large_ext(fq_poly_t beta, slong z_degree) {

	fq_t temp_coeff;
	fq_init(temp_coeff, ctx);

	slong beta_degree = fq_poly_degree(beta, ctx);

	for (slong i = 0; i <= beta_degree; i++) {

		fq_poly_get_coeff(temp_coeff, beta, i, ctx);
		compute_simple_compose(temp_coeff, z_degree);
		fq_poly_set_coeff(beta, i, temp_coeff, ctx);
	}

	fq_clear(temp_coeff, ctx);
}

/**
 * Given $\beta = \sum_i c_i(z)x^i$, computes $\beta = \sum_i c_i(z^{z_degree})x^i$.
 */
void CyclotomicExtRthRoot::compute_beta_coeffs(fq_poly_t beta, slong zeta_degree) {
	Util util;
	if (util.is_small_cyclotomic_ext(r, fq_ctx_prime(ctx)))
		compute_beta_coeffs_small_ext(beta, zeta_degree);
	else
		compute_beta_coeffs_large_ext(beta, zeta_degree);
}
/**
 * Computes $f(z^{z_degree))$.
 */
void CyclotomicExtRthRoot::compute_simple_compose(fq_t f, slong z_degree) {

	fmpz_t temp_coeff1;
	fmpz_t temp_coeff2;
	fmpz_mod_poly_t result;

	fmpz_init(temp_coeff1);
	fmpz_init(temp_coeff2);
	fmpz_mod_poly_init(result, fq_ctx_prime(ctx));

	fmpz_poly_get_coeff_fmpz(temp_coeff1, f, 0);
	fmpz_mod_poly_set_coeff_fmpz(result, 0, temp_coeff1);

	slong degree = fmpz_poly_degree(f);
	slong result_degree = z_degree;

	for (slong i = 1; i <= degree; i++) {

		fmpz_poly_get_coeff_fmpz(temp_coeff1, f, i);
		fmpz_mod_poly_get_coeff_fmpz(temp_coeff2, result, result_degree);

		fmpz_add(temp_coeff1, temp_coeff1, temp_coeff2);
		fmpz_mod_poly_set_coeff_fmpz(result, result_degree, temp_coeff1);

		result_degree = (result_degree + z_degree) % r;
	}

	fmpz_mod_poly_rem(result, result, ctx->modulus);
	fmpz_mod_poly_get_fmpz_poly(f, result);

	fmpz_clear(temp_coeff1);
	fmpz_clear(temp_coeff2);
	fmpz_mod_poly_clear(result);
}

/**
 * Computes $f(x^{x_degree))$.
 */
void CyclotomicExtRthRoot::compute_simple_compose(fq_poly_t f, slong x_degree, const fq_poly_t modulus) {

	fq_t temp_coeff1;
	fq_t temp_coeff2;
	fq_t a_power;
	fq_poly_t result;

	fq_init(temp_coeff1, ctx);
	fq_init(temp_coeff2, ctx);
	fq_init(a_power, ctx);
	fq_poly_init(result, ctx);

	fq_poly_get_coeff(temp_coeff1, f, 0, ctx);
	fq_poly_set_coeff(result, 0, temp_coeff1, ctx);
	fq_one(a_power, ctx);

	slong degree = fq_poly_degree(f, ctx);
	slong result_degree = x_degree;

	for (slong i = 1; i <= degree; i++) {

		fq_poly_get_coeff(temp_coeff1, f, i, ctx);

		if (result_degree >= r) {
			result_degree %= r;
			fq_mul(a_power, a_power, a, ctx);
		}

		fq_mul(temp_coeff1, temp_coeff1, a_power, ctx);
		fq_poly_get_coeff(temp_coeff2, result, result_degree, ctx);
		fq_add(temp_coeff1, temp_coeff1, temp_coeff2, ctx);
		fq_poly_set_coeff(result, result_degree, temp_coeff1, ctx);

		result_degree += x_degree;
	}

	fq_poly_set(f, result, ctx);

	fq_clear(temp_coeff1, ctx);
	fq_clear(temp_coeff2, ctx);
	fq_clear(a_power, ctx);
	fq_poly_clear(result, ctx);
}

/**
 * Compute the initial values $\beta_{init} = \alpha^p$, and $\xi_{init} = x^p$
 */
void CyclotomicExtRthRoot::compute_initials(const fq_poly_t alpha, const fq_poly_t modulus) {

	fmpz_t quotient;
	slong remainder;
	fq_t temp_coeff;

	fmpz_init(quotient);
	fq_init(temp_coeff, ctx);

	// compute xi_init = x^p
	remainder = fmpz_fdiv_ui(fq_ctx_prime(ctx), r);
	fmpz_sub_ui(quotient, fq_ctx_prime(ctx), remainder);
	fmpz_divexact_ui(quotient, quotient, r);
	fq_pow(temp_coeff, a, quotient, ctx);
	fq_poly_zero(xi_init, ctx);
	fq_poly_set_coeff(xi_init, remainder, temp_coeff, ctx);

	// compute beta_init = alpha^p
	fq_poly_set(beta_init, alpha, ctx);
	compute_beta(beta_init, xi_init, remainder, modulus);

	fmpz_clear(quotient);
	fq_clear(temp_coeff, ctx);
}

/**
 * Given a proper factor of y^r - a, computes an rth root of a. Let d = deg(factor), 
 * and let b be constent term of factor. Then (d, r) = 1, so there are integers u, v
 * such that ud + vr = 1. This method computes (-1)^{du}b^u*a^v which is an rth root of a. 
 */
void CyclotomicExtRthRoot::compute_rth_root_from_factor(fq_t root, const fq_poly_t factor) {
	slong d = fq_poly_degree(factor, ctx);
	ulong u = 0;
	ulong v = 0;

	// compute u, v such that ud - vr = 1
	n_xgcd(&u, &v, d, r);

	fq_t temp_coeff;
	fq_t a_inv_pow;
	fq_init(temp_coeff, ctx);
	fq_init(a_inv_pow, ctx);

	fq_poly_get_coeff(temp_coeff, factor, 0, ctx);
	fq_pow_ui(temp_coeff, temp_coeff, u, ctx);
	fq_inv(a_inv_pow, a, ctx);
	fq_pow_ui(a_inv_pow, a_inv_pow, v, ctx);
	fq_mul(a_inv_pow, a_inv_pow, temp_coeff, ctx);

	if ((d % 2) && (u % 2))
		fq_neg(a_inv_pow, a_inv_pow, ctx);

	fq_set(root, a_inv_pow, ctx);

	fq_clear(temp_coeff, ctx);
	fq_clear(a_inv_pow, ctx);
}

void CyclotomicExtRthRoot::compute_rth_root(fq_t root, fq_poly_t f) {

	fq_poly_t alpha;
	fq_poly_t trace;
	fmpz_t power;

	fq_poly_init(alpha, ctx);
	fq_poly_init(trace, ctx);

	flint_rand_t state;
	flint_randinit(state);

	fq_poly_t f_inv_rev;
	fq_poly_init(f_inv_rev, ctx);
	fq_poly_reverse(f_inv_rev, f, fq_poly_length(f, ctx), ctx);
	fq_poly_inv_series_newton(f_inv_rev, f_inv_rev, fq_poly_length(f, ctx), ctx);

	while (true) {

		fq_poly_zero(trace, ctx);

		// find a nonzero trace
		while (fq_poly_is_zero(trace, ctx)) {
			fq_poly_randtest(alpha, state, fq_poly_degree(f, ctx), ctx);
			compute_trace(trace, alpha, f);
		}

		// compute power = (p - 1) / 2
		fmpz_sub_ui(power, fq_ctx_prime(ctx), 1);
		fmpz_divexact_ui(power, power, 2);

		// compute trace = trace^power
		fq_poly_powmod_fmpz_binexp_preinv(trace, trace, power, f, f_inv_rev, ctx);

		fq_poly_one(alpha, ctx);
		fq_poly_sub(trace, trace, alpha, ctx);
		fq_poly_gcd(alpha, trace, f, ctx);

		// check if alpha is a proper factor of f
		if (0 < fq_poly_degree(alpha, ctx) && fq_poly_degree(alpha, ctx) < fq_poly_degree(f, ctx)) {
			compute_rth_root_from_factor(root, alpha);
			break;
		}
	}

	fq_poly_clear(alpha, ctx);
	fq_poly_clear(trace, ctx);
	fq_poly_clear(f_inv_rev, ctx);
	fmpz_clear(power);
	flint_randclear(state);
}

/**
 * Supports aliasing.
 */
void CyclotomicExtRthRoot::compute_rth_root(fq_t root, const fq_t a, slong r, const fq_ctx_t ctx) {

	this->ctx = ctx;
	this->r = r;
	fq_init(this->a, ctx);
	fq_set(this->a, a, ctx);
	fq_poly_init(beta_init, ctx);
	fq_poly_init(xi_init, ctx);

	fq_poly_t f;
	fq_t temp_coeff;

	fq_poly_init(f, ctx);
	fq_init(temp_coeff, ctx);

	fq_one(temp_coeff, ctx);
	fq_poly_set_coeff(f, r, temp_coeff, ctx);
	fq_neg(temp_coeff, a, ctx);
	fq_poly_set_coeff(f, 0, temp_coeff, ctx);

	compute_rth_root(root, f);

	fq_poly_clear(f, ctx);
	fq_clear(temp_coeff, ctx);
	fq_clear(this->a, ctx);
	fq_poly_clear(beta_init, ctx);
	fq_poly_clear(xi_init, ctx);
	this->ctx = NULL;
}
