/*
 * ffIsomorphism.cpp
 *
 *  Created on: Dec 3, 2013
 *      Author: javad
 */

#include "ff_isom.h"
#include "fq_nmod_poly_eval.h"
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
void FFIsomorphism::compute_semi_trace_small_ext(fq_nmod_poly_t delta, fq_nmod_t xi, slong n, const fq_nmod_ctx_t ctx,
		const fq_nmod_poly_t modulus) {

	if (n == 1) {
		fq_nmod_poly_set(delta, delta_init, ctx);
		fq_nmod_set(xi, xi_init, ctx);

		return;
	}

	slong z_degree = 0;
	fq_nmod_poly_t temp_delta;
	fq_nmod_t temp_xi;

	fq_nmod_poly_init(temp_delta, ctx);
	fq_nmod_init(temp_xi, ctx);

	if (n % 2 == 0) {

		compute_semi_trace_small_ext(delta, xi, n / 2, ctx, modulus);
		fq_nmod_poly_set(temp_delta, delta, ctx);
		fq_nmod_set(temp_xi, xi, ctx);
		z_degree = fq_nmod_ctx_degree(ctx) - n / 2;

	} else {

		compute_semi_trace_small_ext(delta, xi, n - 1, ctx, modulus);
		fq_nmod_poly_set(temp_delta, delta_init, ctx);
		fq_nmod_set(temp_xi, xi_init, ctx);
		z_degree = fq_nmod_ctx_degree(ctx) - 1;

	}

	compute_delta(delta, temp_xi, z_degree, ctx, modulus);
	fq_nmod_poly_add(delta, delta, temp_delta, ctx);

	compute_xi(xi, temp_xi, ctx);

	fq_nmod_poly_clear(temp_delta, ctx);
	fq_nmod_clear(temp_xi, ctx);
}

/**
 * Computes the value $\delta = a + z^{r - 1}a^{p^s} + z^{r - 2}a^{p^{2s}} + \cdots + 
 * za^{p^{(r - 1)s}$ where $a \in \mathbb{F}_p[x][z]$, and $r$ is the degree of
 * the extension {@code ctx} of the prime field.
 */
void FFIsomorphism::compute_semi_trace_small_ext(fq_nmod_poly_t theta, const fq_nmod_poly_t a, const fq_nmod_ctx_t ctx,
		const fq_nmod_poly_t modulus) {

	fq_nmod_t xi;
	fq_nmod_init(xi, ctx);

	fq_nmod_poly_set(delta_init, a, ctx);
	compute_semi_trace_small_ext(theta, xi, fq_nmod_ctx_degree(ctx), ctx, modulus);

	fq_nmod_clear(xi, ctx);
}

/**
 * Computes {@code xi} = {@code xi}({@code old_xi}).
 */
void FFIsomorphism::compute_xi(fq_nmod_t xi, const fq_nmod_t old_xi, const fq_nmod_ctx_t ctx) {
	nmod_poly_compose_mod(xi, xi, old_xi, ctx->modulus);
}

/**
 * Given $\delta = \sum_i c_i(x)z^i$, this method computes 
 * $z^{z_degree}\delta(\xi) = z^{z_degree}\sum_i c_i(\xi)z^i$. 
 */
void FFIsomorphism::compute_delta(fq_nmod_poly_t delta, const fq_nmod_t xi, slong z_degree, const fq_nmod_ctx_t ctx,
		const fq_nmod_poly_t modulus) {
	// temporary elements for composition
	fq_nmod_t temp_comp;

	// temporary element for coefficients
	fq_nmod_t temp_coeff;

	fq_nmod_init(temp_comp, ctx);
	fq_nmod_init(temp_coeff, ctx);
	fq_nmod_set(temp_comp, xi, ctx);

	slong delta_degree = fq_nmod_poly_degree(delta, ctx);

	// do modular composition for each coefficient
	for (slong i = 0; i <= delta_degree; i++) {

		fq_nmod_poly_get_coeff(temp_coeff, delta, i, ctx);
		nmod_poly_compose_mod(temp_coeff, temp_coeff, temp_comp, ctx->modulus);
		fq_nmod_poly_set_coeff(delta, i, temp_coeff, ctx);
	}

	// compute z^{z_degree}delta
	// TODO slightly better complexity can be obtained  
	fq_nmod_poly_shift_left(delta, delta, z_degree, ctx);
	fq_nmod_poly_rem(delta, delta, modulus, ctx);

	fq_nmod_clear(temp_comp, ctx);
	fq_nmod_clear(temp_coeff, ctx);
}

/**
 * Computes $\xi_{init} = x^{p^s}$
 */
void FFIsomorphism::compute_xi_init(const fq_nmod_ctx_t ctx, slong s) {
	fq_nmod_t temp_comp1;
	fq_nmod_t temp_comp2;

	fq_nmod_init(temp_comp1, ctx);
	fq_nmod_init(temp_comp2, ctx);

	fq_nmod_zero(temp_comp1, ctx);
	// set temp_comp1 to x
	nmod_poly_set_coeff_ui(temp_comp1, 1, 1);
	// compute x^p
	fq_nmod_pow_ui(temp_comp1, temp_comp1, ctx->modulus->mod.n, ctx);
	fq_nmod_set(temp_comp2, temp_comp1, ctx);

	// reverse the bits of s, needed for binary-powering
	slong bit_length = n_flog(s, 2) + 1;
	slong n = n_revbin(s, bit_length);

	// compute x^{p^s} using a binary-powering scheme
	for (slong i = 1; i < bit_length; i++) {
		n >>= 1;
		nmod_poly_compose_mod(temp_comp2, temp_comp2, temp_comp2, ctx->modulus);
		if (n & 1)
			nmod_poly_compose_mod(temp_comp2, temp_comp2, temp_comp1, ctx->modulus);
	}

	fq_nmod_set(xi_init, temp_comp2, ctx);

	fq_nmod_clear(temp_comp1, ctx);
	fq_nmod_clear(temp_comp2, ctx);
}

/**
 * Computes the value $\theta = \alpha + z^{r - 1}\alpha^{p^s} + z^{r - 2}\alpha^{p^{2s}} + \cdots + 
 * z\alpha^{p^{(r - 1)s}$ where $\alpha \in \mathbb{F}_p[x]$, and $r$ is the degree of
 * the extension {@code ctx} of the prime field.
 */
void FFIsomorphism::compute_semi_trace_large_ext(fq_nmod_poly_t theta, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx,
		const fq_nmod_poly_t modulus) {

	slong degree = fq_nmod_ctx_degree(ctx);

	fq_nmod_t *frobenius = new fq_nmod_t[degree];
	for (slong i = 0; i < degree; i++)
		fq_nmod_init(frobenius[i], ctx);

	iterated_frobenius(frobenius, alpha, ctx, fq_nmod_poly_degree(modulus, ctx));

	fq_nmod_poly_zero(theta, ctx);
	fq_nmod_poly_set_coeff(theta, 0, frobenius[0], ctx);
	for (slong i = 1; i < degree; i++) {
		fq_nmod_poly_set_coeff(theta, i, frobenius[degree - i], ctx);
	}

	fq_nmod_poly_rem(theta, theta, modulus, ctx);

	for (slong i = 0; i < degree; i++)
		fq_nmod_clear(frobenius[i], ctx);
	delete[] frobenius;
}

/**
 * Given an elemenet {@code alpha} in the field {@code ctx}, this method computes
 * $\alpha^{p^{is}}$ for all $i = 0, \dots, r - 1$. The algorithm used is Algorithm 3.1
 * form Von Zur Gathen and Shoup, 1992. More precisely, this method implements the
 * case {q = p, R = ctx, t = p^s, m = r - 1} of the algorithm in the paper.
 */
void FFIsomorphism::iterated_frobenius(fq_nmod_t *result, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx, slong s) {
	fq_nmodPolyEval fq_nmodPolyEval;

	fq_nmod_poly_t temp;
	fq_nmod_poly_init(temp, ctx);

	fq_nmod_zero(result[0], ctx);
	// set result[0] to x
	nmod_poly_set_coeff_ui(result[0], 1, 1);
	fq_nmod_set(result[1], xi_init, ctx);

	slong l = n_clog(fq_nmod_ctx_degree(ctx) - 1, 2);
	slong base = 0;
	slong length = 0;

	for (slong i = 1; i <= l; i++) {

		base = 1 << (i - 1);

		// build the polynomial for multipoint evaluation
		convert(temp, result[base], ctx);
		
		// make sure we stay in the bound
		if (2 * base < fq_nmod_ctx_degree(ctx))
			length = base;
		else
			length = fq_nmod_ctx_degree(ctx) - base - 1;

		fq_nmodPolyEval.multipoint_eval(result + base + 1, temp, result + 1, length, ctx);
	}

	// check the trivial case of alpha = x
	if (!fq_nmod_equal(result[0], alpha, ctx)) {
		
		// build the polynomial for multipoint evaluation of alpha
		convert(temp, alpha, ctx);

		fq_nmodPolyEval.multipoint_eval(result, temp, result, fq_nmod_ctx_degree(ctx), ctx);
	}

	fq_nmod_poly_clear(temp, ctx);
}

/**
 * Buidls a cyclotomic extension of the prime field by computing an irreducible
 * factor of the $r$-th cyclotomic polynomial.
 * 
 * @param modulus	fq_nmod_poly version of the resulting cyclotomic factor
 * @param ctx		the resulting cyclotomic extension
 */
void FFIsomorphism::build_cyclotomic_extension(fq_nmod_poly_t modulus, fq_nmod_ctx_t cyclotomic_ctx) {
	fq_nmod_t cyclo_poly;
	fq_nmod_init(cyclo_poly, ctx_1);

	Util util;
	// degree of the given extension
	slong r = fq_nmod_ctx_degree(ctx_1);
	// degree of the cyclotomic extension
	slong s = util.compute_multiplicative_order(ctx_1->modulus->mod.n, r);

	// build the r-th cyclotomic polynomial
	for (slong i = 0; i < r; i++)
		nmod_poly_set_coeff_ui(cyclo_poly, i, 1);

	// obtain an irreducible factor
	nmod_poly_factor_t factors;
	nmod_poly_factor_init(factors);
	nmod_poly_factor_equal_deg(factors, cyclo_poly, s);
	fq_nmod_ctx_init_modulus(cyclotomic_ctx, &factors->p[0], "z");

	// build the modulus
	convert(modulus, cyclotomic_ctx->modulus, ctx_1);

	nmod_poly_factor_clear(factors);
	nmod_poly_clear(cyclo_poly);
}

/**
 * Coerces the element {@code value} to an element of the field {@code ctx},
 * and store it in {@code result}. It is assumed that the coefficients of {@code value}
 * are in the base field.
 */
void FFIsomorphism::convert(fq_nmod_t result, const fq_nmod_poly_t value, const fq_nmod_ctx_t ctx) {
	fq_nmod_t temp_coeff1;
	mp_limb_t temp_coeff2 = 0;

	fq_nmod_init(temp_coeff1, ctx);
	fq_nmod_zero(result, ctx);

	for (slong i = 0; i <= fq_nmod_poly_degree(value, ctx); i++) {
		fq_nmod_poly_get_coeff(temp_coeff1, value, i, ctx);
		temp_coeff2 = nmod_poly_get_coeff_ui(temp_coeff1, 0);
		nmod_poly_set_coeff_ui(result, i, temp_coeff2);
	}

	fq_nmod_clear(temp_coeff1, ctx);
}

/**
 * Converts the element {@code value} to a polynomial over the field {@code ctx},
 * and store it in {@code result}.
 */
void FFIsomorphism::convert(fq_nmod_poly_t result, const fq_nmod_t value, const fq_nmod_ctx_t ctx) {
	mp_limb_t temp_coeff = 0;
	
	fq_nmod_t temp;
	fq_nmod_init(temp, ctx);

	fq_nmod_poly_zero(result, ctx);

	for (slong i = 0; i <= nmod_poly_degree(value); i++) {
		temp_coeff = nmod_poly_get_coeff_ui(value, i);
		nmod_poly_set_coeff_ui(temp, 0, temp_coeff);
		fq_nmod_poly_set_coeff(result, i, temp, ctx);
	}

	fq_nmod_clear(temp, ctx);
}

/**
 * Compute the isomorphism between the two extensions of the form $\mathbb{F}_p[z][x] / (x^r - \eta)$.
 * The isomorphism is of the form $x \mapsto cx$ for some $c \in \mathbb{F}_p[z]$. 
 */
void FFIsomorphism::compute_middle_isomorphism(fq_nmod_t c, const fq_nmod_poly_t theta_a, const fq_nmod_poly_t theta_b,
		const fq_nmod_poly_t modulus, const fq_nmod_ctx_t cyclotomic_ctx) {

	fq_nmod_poly_t eta;
	fq_nmod_poly_t modulus_inv_rev;

	fq_nmod_poly_init(modulus_inv_rev, ctx_1);
	fq_nmod_poly_init(eta, ctx_1);

	fq_nmod_t temp;
	fq_nmod_init(temp, cyclotomic_ctx);

	fq_nmod_poly_reverse(modulus_inv_rev, modulus, fq_nmod_poly_length(modulus, ctx_1), ctx_1);
	fq_nmod_poly_inv_series_newton(modulus_inv_rev, modulus_inv_rev, fq_nmod_poly_length(modulus, ctx_1), ctx_1);

	// compute theta_a^r
	fq_nmod_poly_powmod_ui_binexp_preinv(eta, theta_a, fq_nmod_ctx_degree(ctx_1), modulus, modulus_inv_rev, ctx_1);
	// now eta is an element in the cyclotomic field
	convert(c, eta, cyclotomic_ctx);

	// compute theta_b^r
	fq_nmod_poly_powmod_ui_binexp_preinv(eta, theta_b, fq_nmod_ctx_degree(ctx_2), modulus, modulus_inv_rev, ctx_2);
	convert(temp, eta, cyclotomic_ctx);

	// compute c = theta_a^r / theta_b^r
	fq_nmod_inv(temp, temp, cyclotomic_ctx);
	fq_nmod_mul(c, c, temp, cyclotomic_ctx);

	// compute c^{1 / r}
	CyclotomicExtRthRoot cyclotomicExtRthRoot;
	cyclotomicExtRthRoot.compute_rth_root(c, c, fq_nmod_ctx_degree(ctx_1), cyclotomic_ctx);

	fq_nmod_poly_clear(eta, ctx_1);
	fq_nmod_poly_clear(modulus_inv_rev, ctx_1);
	fq_nmod_clear(temp, cyclotomic_ctx);
}

/**
 * Computes a semi-trace. This methods decides the proper semi-trace computation approach
 * based on the degree of the auxiliary cyclotomic extension. 
 * 
 * @param theta		the resulting semi-trace
 */
void FFIsomorphism::compute_semi_trace(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus) {

	flint_rand_t state;
	flint_randinit(state);
	fq_nmod_poly_zero(theta, ctx);

	Util util;
	slong degree = fq_nmod_ctx_degree(ctx);
	slong s = fq_nmod_poly_degree(modulus, ctx);

	// computing xi_init = x^{p^s}
	compute_xi_init(ctx, s);

	if (util.is_small_cyclotomic_ext(degree, ctx->modulus->mod.n)) {
		fq_nmod_poly_t alpha;
		fq_nmod_poly_init(alpha, ctx);

		while (fq_nmod_poly_is_zero(theta, ctx)) {
			fq_nmod_poly_randtest_not_zero(alpha, state, s, ctx);
			compute_semi_trace_small_ext(theta, alpha, ctx, modulus);
		}

		fq_nmod_poly_clear(alpha, ctx);
		flint_randclear(state);

		return;
	}

	fq_nmod_t alpha;
	fq_nmod_init(alpha, ctx);

	// try alpha = x first
	nmod_poly_set_coeff_ui(alpha, 1, 1);
	compute_semi_trace_large_ext(theta, alpha, ctx, modulus);
	
	// if the semi trace of x is zero we try random cases
	while (fq_nmod_poly_is_zero(theta, ctx)) {
		fq_nmod_randtest_not_zero(alpha, state, ctx);
		compute_semi_trace_large_ext(theta, alpha, ctx, modulus);
	}

	fq_nmod_clear(alpha, ctx);
	flint_randclear(state);
}

/**
 * Compute the isomorphism between the cyclotomic extensions of {@code ctx_1}, {@code ctx_2}.
 * The resulting isomorphism is $f \mapsto f_{image}$.
 */
void FFIsomorphism::compute_extension_isomorphism(fq_nmod_poly_t f, fq_nmod_poly_t f_image) {
	fq_nmod_ctx_t cyclotomic_ctx;
	fq_nmod_poly_t modulus;
	fq_nmod_poly_init(modulus, ctx_1);
	build_cyclotomic_extension(modulus, cyclotomic_ctx);

	compute_semi_trace(f, ctx_1, modulus);
	compute_semi_trace(f_image, ctx_2, modulus);

	fq_nmod_t c;
	fq_nmod_poly_t c_temp;
	fq_nmod_init(c, cyclotomic_ctx);
	fq_nmod_poly_init(c_temp, ctx_2);

	compute_middle_isomorphism(c, f, f_image, modulus, cyclotomic_ctx);

	convert(c_temp, c, ctx_2);
	fq_nmod_poly_mulmod(f_image, f_image, c_temp, modulus, ctx_2);

	fq_nmod_poly_clear(modulus, ctx_1);
	fq_nmod_poly_clear(c_temp, ctx_2);
	fq_nmod_clear(c, cyclotomic_ctx);
	fq_nmod_ctx_clear(cyclotomic_ctx);
}

/**
 * Builds an isomorphism
 * h: ctx_1 --> ctx_2
 *        x --> x_image
 */
void FFIsomorphism::build_isomorphism() {
	fq_nmod_poly_init(delta_init, ctx_1);
	fq_nmod_init(xi_init, ctx_1);

	fq_nmod_poly_t f;
	fq_nmod_poly_t f_image;

	fq_nmod_poly_init(f, ctx_1);
	fq_nmod_poly_init(f_image, ctx_2);

	compute_extension_isomorphism(f, f_image);

	fq_nmod_t f0;
	fq_nmod_t h0;

	fq_nmod_init(f0, ctx_1);
	fq_nmod_init(h0, ctx_2);

	fq_nmod_poly_get_coeff(f0, f, 0, ctx_1);
	fq_nmod_poly_get_coeff(h0, f_image, 0, ctx_2);

	// now we have an isomorphism between 
	// h: ctx_1 --> ctx_2
	//        f0 --> h0
	// and we want to compute an isomorphism
	// h: ctx_1 --> ctx_2
	//        x --> g
	// for some g

	fq_nmod_t x;
	fq_nmod_init(x, ctx_1);
	nmod_poly_set_coeff_ui(x, 1, 1);
	
	FFIsomBaseChange ffIsomBaseChange;
	ffIsomBaseChange.change_basis(x_image, f0, x, ctx_1->modulus);
	nmod_poly_compose_mod(x_image, x_image, h0, ctx_2->modulus);

	fq_nmod_poly_clear(f, ctx_1);
	fq_nmod_poly_clear(f_image, ctx_2);
	fq_nmod_clear(f0, ctx_1);
	fq_nmod_clear(h0, ctx_2);
	fq_nmod_clear(x, ctx_2);
	fq_nmod_poly_clear(delta_init, ctx_1);
	fq_nmod_clear(xi_init, ctx_1);
}

void FFIsomorphism::compute_image(nmod_poly_t image, const nmod_poly_t f) {
	nmod_poly_compose_mod(image, f, x_image, ctx_2->modulus);
}

FFIsomorphism::FFIsomorphism(const nmod_poly_t modulus1, const nmod_poly_t modulus2) {
	nmod_poly_t tempf1;
	nmod_poly_t tempf2;
	nmod_poly_init(tempf1, modulus1->mod.n);
	nmod_poly_init(tempf2, modulus2->mod.n);
	nmod_poly_set(tempf1, modulus1);
	nmod_poly_set(tempf2, modulus2);

	fq_nmod_ctx_init_modulus(ctx_1, tempf1, "x");
	fq_nmod_ctx_init_modulus(ctx_2, tempf2, "x");
	nmod_poly_init(x_image, modulus2->mod.n);

	build_isomorphism();

	nmod_poly_clear(tempf1);
	nmod_poly_clear(tempf2);
}

FFIsomorphism::~FFIsomorphism() {
	fq_nmod_ctx_clear(ctx_1);
	fq_nmod_ctx_clear(ctx_2);
	nmod_poly_clear(x_image);
}
