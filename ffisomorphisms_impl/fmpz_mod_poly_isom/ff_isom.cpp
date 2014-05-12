/*
 * ffIsomorphism.cpp
 *
 *  Created on: Dec 3, 2013
 *      Author: javad
 */

#include "ff_isom.h"
#include "fq_poly_eval.h"
#include "cyclotomic_ext_rth_root.h"
#include "ff_isom_base_change.h"
#include "util.h"
#include <iostream>
#include <flint/profiler.h>

using namespace std;

/**
 * Computes the value $\delta_n = a + z^{r - 1}a^{p^s} + z^{r - 2}a^{p^{2s}} + \cdots + 
 * z^{r - n + 1}a^{p^{(n - 1)s}$ where $a \in \mathbb{F}_p[x][z]$ is such that $\delta_init = a$.
 * After recursion level i, $\delta_n = a + z^{r - 1}a^{p^s} + z^{r - 2}a^{p^{2s}} + \cdots + 
 * z^{r - i + 1}a^{p^{(i - 1)s}$ and $\xi_i = x^{p^{is}}$.
 */
void FFIsomorphism::compute_semi_trace_small_ext(fq_poly_t delta, fq_t xi, slong n, const fq_ctx_t ctx,
		const fq_poly_t modulus) {

	if (n == 1) {
		fq_poly_set(delta, delta_init, ctx);
		fq_set(xi, xi_init, ctx);

		return;
	}

	slong z_degree = 0;
	fq_poly_t temp_delta;
	fq_t temp_xi;

	fq_poly_init(temp_delta, ctx);
	fq_init(temp_xi, ctx);

	if (n % 2 == 0) {

		compute_semi_trace_small_ext(delta, xi, n / 2, ctx, modulus);
		fq_poly_set(temp_delta, delta, ctx);
		fq_set(temp_xi, xi, ctx);
		z_degree = fq_ctx_degree(ctx) - n / 2;

	} else {

		compute_semi_trace_small_ext(delta, xi, n - 1, ctx, modulus);
		fq_poly_set(temp_delta, delta_init, ctx);
		fq_set(temp_xi, xi_init, ctx);
		z_degree = fq_ctx_degree(ctx) - 1;

	}

	compute_delta(delta, temp_xi, z_degree, ctx, modulus);
	fq_poly_add(delta, delta, temp_delta, ctx);

	compute_xi(xi, temp_xi, ctx);

	fq_poly_clear(temp_delta, ctx);
	fq_clear(temp_xi, ctx);
}

/**
 * Computes the value $\delta = a + z^{r - 1}a^{p^s} + z^{r - 2}a^{p^{2s}} + \cdots + 
 * za^{p^{(r - 1)s}$ where $a \in \mathbb{F}_p[x][z]$, and $r$ is the degree of
 * the extension {@code ctx} of the prime field.
 */
void FFIsomorphism::compute_semi_trace_small_ext(fq_poly_t theta, const fq_poly_t a, const fq_ctx_t ctx,
		const fq_poly_t modulus) {

	fq_t xi;
	fq_init(xi, ctx);

	fq_poly_set(delta_init, a, ctx);
	compute_semi_trace_small_ext(theta, xi, fq_ctx_degree(ctx), ctx, modulus);

	fq_clear(xi, ctx);
}

/**
 * Computes {@code xi} = {@code xi}({@code old_xi}).
 */
void FFIsomorphism::compute_xi(fq_t xi, const fq_t old_xi, const fq_ctx_t ctx) {
	// temporary elements for composition
	fmpz_mod_poly_t temp_comp1;
	fmpz_mod_poly_t temp_comp2;

	fmpz_mod_poly_init(temp_comp1, fq_ctx_prime(ctx));
	fmpz_mod_poly_init(temp_comp2, fq_ctx_prime(ctx));

	fmpz_mod_poly_set_fmpz_poly(temp_comp1, xi);
	fmpz_mod_poly_set_fmpz_poly(temp_comp2, old_xi);

	fmpz_mod_poly_compose_mod(temp_comp1, temp_comp1, temp_comp2, ctx->modulus);
	fmpz_mod_poly_get_fmpz_poly(xi, temp_comp1);

	fmpz_mod_poly_clear(temp_comp1);
	fmpz_mod_poly_clear(temp_comp2);
}

/**
 * Given $\delta = \sum_i c_i(x)z^i$, this method computes 
 * $z^{z_degree}\delta(\xi) = z^{z_degree}\sum_i c_i(\xi)z^i$. 
 */
void FFIsomorphism::compute_delta(fq_poly_t delta, const fq_t xi, slong z_degree, const fq_ctx_t ctx,
		const fq_poly_t modulus) {
	// temporary elements for composition
	fmpz_mod_poly_t temp_comp1;
	fmpz_mod_poly_t temp_comp2;

	// temporary element for coefficients
	fq_t temp_coeff;

	fmpz_mod_poly_init(temp_comp1, fq_ctx_prime(ctx));
	fmpz_mod_poly_init(temp_comp2, fq_ctx_prime(ctx));
	fq_init(temp_coeff, ctx);
	fmpz_mod_poly_set_fmpz_poly(temp_comp1, xi);

	slong delta_degree = fq_poly_degree(delta, ctx);

	// do modular composition for each coefficient
	for (slong i = 0; i <= delta_degree; i++) {

		fq_poly_get_coeff(temp_coeff, delta, i, ctx);
		fmpz_mod_poly_set_fmpz_poly(temp_comp2, temp_coeff);
		fmpz_mod_poly_compose_mod(temp_comp2, temp_comp2, temp_comp1, ctx->modulus);
		fmpz_mod_poly_get_fmpz_poly(temp_coeff, temp_comp2);
		fq_poly_set_coeff(delta, i, temp_coeff, ctx);
	}

	// compute z^{z_degree}delta
	// TODO slightly better complexity can be obtained  
	fq_poly_shift_left(delta, delta, z_degree, ctx);
	fq_poly_rem(delta, delta, modulus, ctx);

	fmpz_mod_poly_clear(temp_comp1);
	fmpz_mod_poly_clear(temp_comp2);
	fq_clear(temp_coeff, ctx);
}

/**
 * Computes $\xi_{init} = x^{p^s}$
 */
void FFIsomorphism::compute_xi_init(const fq_ctx_t ctx, slong s) {
	fmpz_mod_poly_t temp_comp1;
	fmpz_mod_poly_t temp_comp2;

	fmpz_mod_poly_init(temp_comp1, fq_ctx_prime(ctx));
	fmpz_mod_poly_init(temp_comp2, fq_ctx_prime(ctx));

	fmpz_mod_poly_zero(temp_comp1);
	// set temp_comp1 to x
	fmpz_mod_poly_set_coeff_ui(temp_comp1, 1, 1);
	// compute x^p
	fmpz_mod_poly_powmod_fmpz_binexp(temp_comp1, temp_comp1, fq_ctx_prime(ctx), ctx->modulus);
	fmpz_mod_poly_set(temp_comp2, temp_comp1);

	// reverse the bits of s, needed for binary-powering
	slong bit_length = n_flog(s, 2) + 1;
	slong n = n_revbin(s, bit_length);

	// compute x^{p^s} using a binary-powering scheme
	for (slong i = 1; i < bit_length; i++) {
		n >>= 1;
		fmpz_mod_poly_compose_mod(temp_comp2, temp_comp2, temp_comp2, ctx->modulus);
		if (n & 1)
			fmpz_mod_poly_compose_mod(temp_comp2, temp_comp2, temp_comp1, ctx->modulus);
	}

	fmpz_mod_poly_get_fmpz_poly(xi_init, temp_comp2);

	fmpz_mod_poly_clear(temp_comp1);
	fmpz_mod_poly_clear(temp_comp2);
}

/**
 * Computes the value $\theta = \alpha + z^{r - 1}\alpha^{p^s} + z^{r - 2}\alpha^{p^{2s}} + \cdots + 
 * z\alpha^{p^{(r - 1)s}$ where $\alpha \in \mathbb{F}_p[x]$, and $r$ is the degree of
 * the extension {@code ctx} of the prime field.
 */
void FFIsomorphism::compute_semi_trace_large_ext(fq_poly_t theta, const fq_t alpha, const fq_ctx_t ctx,
		const fq_poly_t modulus) {

	slong degree = fq_ctx_degree(ctx);

	fq_t *frobenius = new fq_t[degree];
	for (slong i = 0; i < degree; i++)
		fq_init(frobenius[i], ctx);

	iterated_frobenius(frobenius, alpha, ctx, fq_poly_degree(modulus, ctx));

	fq_poly_zero(theta, ctx);
	fq_poly_set_coeff(theta, 0, frobenius[0], ctx);
	for (slong i = 1; i < degree; i++) {
		fq_poly_set_coeff(theta, i, frobenius[degree - i], ctx);
	}

	fq_poly_rem(theta, theta, modulus, ctx);

	for (slong i = 0; i < degree; i++)
		fq_clear(frobenius[i], ctx);
	delete[] frobenius;
}

/**
 * Given an elemenet {@code alpha} in the field {@code ctx}, this method computes
 * $\alpha^{p^{is}}$ for all $i = 0, \dots, r - 1$. The algorithm used is Algorithm 3.1
 * form Von Zur Gathen and Shoup, 1992. More precisely, this method implements the
 * case {q = p, R = ctx, t = p^s, m = r - 1} of the algorithm in the paper.
 */
void FFIsomorphism::iterated_frobenius(fq_t *result, const fq_t alpha, const fq_ctx_t ctx, slong s) {
	FqPolyEval fqPolyEval;

	fq_poly_t temp;
	fq_poly_init(temp, ctx);

	fq_zero(result[0], ctx);
	// set result[0] to x
	fmpz_poly_set_coeff_ui(result[0], 1, 1);
	fq_set(result[1], xi_init, ctx);

	slong l = n_clog(fq_ctx_degree(ctx) - 1, 2);
	slong base = 0;
	slong length = 0;

	for (slong i = 1; i <= l; i++) {

		base = 1 << (i - 1);

		// build the polynomial for multipoint evaluation
		convert(temp, result[base], ctx);

		// make sure we stay in the bound
		if (2 * base < fq_ctx_degree(ctx))
			length = base;
		else
			length = fq_ctx_degree(ctx) - base - 1;

		fqPolyEval.multipoint_eval(result + base + 1, temp, result + 1, length, ctx);
	}

	// check the trivial case of alpha = x
	if (!fq_equal(result[0], alpha, ctx)) {
		
		// build the polynomial for multipoint evaluation of alpha
		convert(temp, alpha, ctx);

		fqPolyEval.multipoint_eval(result, temp, result, fq_ctx_degree(ctx), ctx);
	}

	fq_poly_clear(temp, ctx);
}

/**
 * Buidls a cyclotomic extension of the prime field by computing an irreducible
 * factor of the $r$-th cyclotomic polynomial.
 * 
 * @param modulus	fq_poly version of the resulting cyclotomic factor
 * @param ctx		the resulting cyclotomic extension
 */
void FFIsomorphism::build_cyclotomic_extension(fq_poly_t modulus, fq_ctx_t cyclotomic_ctx) {
	fmpz_mod_poly_t cyclo_poly;
	fmpz_mod_poly_init(cyclo_poly, fq_ctx_prime(ctx_1));

	Util util;
	// degree of the given extension
	slong r = fq_ctx_degree(ctx_1);
	// degree of the cyclotomic extension
	slong s = util.compute_multiplicative_order(fq_ctx_prime(ctx_1), r);

	// build the r-th cyclotomic polynomial
	for (slong i = 0; i < r; i++)
		fmpz_mod_poly_set_coeff_ui(cyclo_poly, i, 1);

	// obtain an irreducible factor
	fmpz_mod_poly_factor_t factors;
	fmpz_mod_poly_factor_init(factors);
	fmpz_mod_poly_factor_equal_deg(factors, cyclo_poly, s);
	fq_ctx_init_modulus(cyclotomic_ctx, &factors->poly[0], "z");

	fmpz_t temp_coeff;
	fmpz_init(temp_coeff);
	
	// build the modulus
	for (slong i = 0; i <= s; i++) {
		fmpz_mod_poly_get_coeff_fmpz(temp_coeff, cyclotomic_ctx->modulus, i);
		fq_poly_set_coeff_fmpz(modulus, i, temp_coeff, ctx_1);
	}

	fmpz_mod_poly_factor_clear(factors);
	fmpz_mod_poly_clear(cyclo_poly);
	fmpz_clear(temp_coeff);
}

/**
 * Coerces the element {@code value} to an element of the field {@code ctx},
 * and store it in {@code result}. It is assumed that the coefficients of {@code value}
 * are in the base field.
 */
void FFIsomorphism::convert(fq_t result, const fq_poly_t value, const fq_ctx_t ctx) {
	fq_t temp_coeff1;
	fmpz_t temp_coeff2;

	fq_init(temp_coeff1, ctx);
	fmpz_init(temp_coeff2);
	fq_zero(result, ctx);

	for (slong i = 0; i <= fq_poly_degree(value, ctx); i++) {
		fq_poly_get_coeff(temp_coeff1, value, i, ctx);
		fmpz_poly_get_coeff_fmpz(temp_coeff2, temp_coeff1, 0);
		fmpz_poly_set_coeff_fmpz(result, i, temp_coeff2);
	}

	fq_clear(temp_coeff1, ctx);
	fmpz_clear(temp_coeff2);
}

/**
 * Converts the element {@code value} to a polynomial over the field {@code ctx},
 * and store it in {@code result}.
 */
void FFIsomorphism::convert(fq_poly_t result, const fq_t value, const fq_ctx_t ctx) {
	fmpz_t temp_coeff;
	fmpz_init(temp_coeff);

	fq_poly_zero(result, ctx);

	for (slong i = 0; i <= fmpz_poly_degree(value); i++) {
		fmpz_poly_get_coeff_fmpz(temp_coeff, value, i);
		fq_poly_set_coeff_fmpz(result, i, temp_coeff, ctx);
	}

	fmpz_clear(temp_coeff);
}

/**
 * Compute the isomorphism between the two extensions of the form $\mathbb{F}_p[z][x] / (x^r - \eta)$.
 * The isomorphism is of the form $x \mapsto cx$ for some $c \in \mathbb{F}_p[z]$. 
 */
void FFIsomorphism::compute_middle_isomorphism(fq_t c, const fq_poly_t theta_a, const fq_poly_t theta_b,
		const fq_poly_t modulus, const fq_ctx_t cyclotomic_ctx) {

	fq_poly_t eta;
	fq_poly_t modulus_inv_rev;

	fq_poly_init(modulus_inv_rev, ctx_1);
	fq_poly_init(eta, ctx_1);

	fq_t temp;
	fq_init(temp, cyclotomic_ctx);

	fq_poly_reverse(modulus_inv_rev, modulus, fq_poly_length(modulus, ctx_1), ctx_1);
	fq_poly_inv_series_newton(modulus_inv_rev, modulus_inv_rev, fq_poly_length(modulus, ctx_1), ctx_1);

	// compute theta_a^r
	fq_poly_powmod_ui_binexp_preinv(eta, theta_a, fq_ctx_degree(ctx_1), modulus, modulus_inv_rev, ctx_1);
	// now eta is an element in the cyclotomic field
	convert(c, eta, cyclotomic_ctx);

	// compute theta_b^r
	fq_poly_powmod_ui_binexp_preinv(eta, theta_b, fq_ctx_degree(ctx_2), modulus, modulus_inv_rev, ctx_2);
	convert(temp, eta, cyclotomic_ctx);

	// compute c = theta_a^r / theta_b^r
	fq_inv(temp, temp, cyclotomic_ctx);
	fq_mul(c, c, temp, cyclotomic_ctx);

	// compute c^{1 / r}
	CyclotomicExtRthRoot cyclotomicExtRthRoot;
	cyclotomicExtRthRoot.compute_rth_root(c, c, fq_ctx_degree(ctx_1), cyclotomic_ctx);

	fq_poly_clear(eta, ctx_1);
	fq_poly_clear(modulus_inv_rev, ctx_1);
	fq_clear(temp, cyclotomic_ctx);
}

/**
 * Computes a semi-trace. This methods decides the proper semi-trace computation approach
 * based on the degree of the auxiliary cyclotomic extension. 
 * 
 * @param theta		the resulting semi-trace
 */
void FFIsomorphism::compute_semi_trace(fq_poly_t theta, const fq_ctx_t ctx, const fq_poly_t modulus) {

	flint_rand_t state;
	flint_randinit(state);
	fq_poly_zero(theta, ctx);

	Util util;
	slong degree = fq_ctx_degree(ctx);
	slong s = fq_poly_degree(modulus, ctx);

	// computing xi_init = x^{p^s}
	compute_xi_init(ctx, s);

	if (util.is_small_cyclotomic_ext(degree, fq_ctx_prime(ctx))) {
		fq_poly_t alpha;
		fq_poly_init(alpha, ctx);

		while (fq_poly_is_zero(theta, ctx)) {
			fq_poly_randtest_not_zero(alpha, state, s, ctx);
			compute_semi_trace_small_ext(theta, alpha, ctx, modulus);
		}

		fq_poly_clear(alpha, ctx);
		flint_randclear(state);

		return;
	}

	fq_t alpha;
	fq_init(alpha, ctx);

	// try alpha = x first
	fmpz_poly_set_coeff_ui(alpha, 1, 1);
	compute_semi_trace_large_ext(theta, alpha, ctx, modulus);
	
	// if the semi trace of x is zero we try random cases
	while (fq_poly_is_zero(theta, ctx)) {
		fq_randtest_not_zero(alpha, state, ctx);
		compute_semi_trace_large_ext(theta, alpha, ctx, modulus);
	}

	fq_clear(alpha, ctx);
	flint_randclear(state);
}

/**
 * Compute the isomorphism between the cyclotomic extensions of {@code ctx_1}, {@code ctx_2}.
 * The resulting isomorphism is $f \mapsto f_{image}$.
 */
void FFIsomorphism::compute_extension_isomorphism(fq_poly_t f, fq_poly_t f_image) {
	fq_ctx_t cyclotomic_ctx;
	fq_poly_t modulus;
	fq_poly_init(modulus, ctx_1);
	build_cyclotomic_extension(modulus, cyclotomic_ctx);

	compute_semi_trace(f, ctx_1, modulus);
	compute_semi_trace(f_image, ctx_2, modulus);

	fq_t c;
	fq_poly_t c_temp;
	fq_init(c, cyclotomic_ctx);
	fq_poly_init(c_temp, ctx_2);

	compute_middle_isomorphism(c, f, f_image, modulus, cyclotomic_ctx);

	convert(c_temp, c, ctx_2);
	fq_poly_mulmod(f_image, f_image, c_temp, modulus, ctx_2);

	fq_poly_clear(modulus, ctx_1);
	fq_poly_clear(c_temp, ctx_2);
	fq_clear(c, cyclotomic_ctx);
	fq_ctx_clear(cyclotomic_ctx);
}

/**
 * Builds an isomorphism
 * h: ctx_1 --> ctx_2
 *        x --> x_image
 */
void FFIsomorphism::build_isomorphism() {
	fq_poly_init(delta_init, ctx_1);
	fq_init(xi_init, ctx_1);

	fq_poly_t f;
	fq_poly_t f_image;

	fq_poly_init(f, ctx_1);
	fq_poly_init(f_image, ctx_2);

	compute_extension_isomorphism(f, f_image);

	fmpz_poly_t temp;
	fmpz_mod_poly_t f0;
	fmpz_mod_poly_t h0;

	fmpz_poly_init(temp);
	fmpz_mod_poly_init(f0, fq_ctx_prime(ctx_1));
	fmpz_mod_poly_init(h0, fq_ctx_prime(ctx_2));

	fq_poly_get_coeff(temp, f, 0, ctx_1);
	fmpz_mod_poly_set_fmpz_poly(f0, temp);
	fq_poly_get_coeff(temp, f_image, 0, ctx_2);
	fmpz_mod_poly_set_fmpz_poly(h0, temp);

	// now we have an isomorphism between 
	// h: ctx_1 --> ctx_2
	//        f0 --> h0
	// and we want to compute an isomorphism
	// h: ctx_1 --> ctx_2
	//        x --> g
	// for some g

	fmpz_mod_poly_t x;
	fmpz_mod_poly_init(x, fq_ctx_prime(ctx_1));
	fmpz_mod_poly_set_coeff_ui(x, 1, 1);
	
	FFIsomBaseChange ffIsomBaseChange;
	ffIsomBaseChange.change_basis(x_image, f0, x, ctx_1->modulus);
	fmpz_mod_poly_compose_mod(x_image, x_image, h0, ctx_2->modulus);

	fq_poly_clear(f, ctx_1);
	fq_poly_clear(f_image, ctx_2);
	fmpz_poly_clear(temp);
	fmpz_mod_poly_clear(f0);
	fmpz_mod_poly_clear(h0);
	fmpz_mod_poly_clear(x);
	fq_poly_clear(delta_init, ctx_1);
	fq_clear(xi_init, ctx_1);
}

void FFIsomorphism::compute_image(fmpz_mod_poly_t image, const fmpz_mod_poly_t f) {
	fmpz_mod_poly_compose_mod(image, f, x_image, ctx_2->modulus);
}

FFIsomorphism::FFIsomorphism(const fmpz_mod_poly_t modulus1, const fmpz_mod_poly_t modulus2) {
	fmpz_mod_poly_t tempf1;
	fmpz_mod_poly_t tempf2;
	fmpz_mod_poly_init(tempf1, &modulus1->p);
	fmpz_mod_poly_init(tempf2, &modulus2->p);
	fmpz_mod_poly_set(tempf1, modulus1);
	fmpz_mod_poly_set(tempf2, modulus2);

	fq_ctx_init_modulus(ctx_1, tempf1, "x");
	fq_ctx_init_modulus(ctx_2, tempf2, "x");
	fmpz_mod_poly_init(x_image, &modulus2->p);

	build_isomorphism();

	fmpz_mod_poly_clear(tempf1);
	fmpz_mod_poly_clear(tempf2);
}

FFIsomorphism::~FFIsomorphism() {
	fq_ctx_clear(ctx_1);
	fq_ctx_clear(ctx_2);
	fmpz_mod_poly_clear(x_image);
}
