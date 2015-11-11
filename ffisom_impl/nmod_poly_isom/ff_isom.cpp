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
#include "nmod_cyclotomic_poly.h"
#include "util.h"
#include <iostream>
#include <flint/profiler.h>

using namespace std;

/**
 * Computes the value $\delta_n = a + z^{r - 1}\sigma(a) + z^{r - 2}\sigma^2(a) + \cdots + 
 * z^{r - n + 1}\sigma^{n - 1}(a)$ where $a \in \mathbb{F}_p[x][z]$ is such that $\delta_init = a$.
 * After recursion level i, $\delta_n = a + z^{r - 1}\sigma(a) + z^{r - 2}\sigma^2(a) + \cdots + 
 * z^{r - i + 1}\sigma^{i - 1}(a)$ and $\xi_i = x^{p^i}$.
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

void FFIsomorphism::compute_semi_trace_trivial_ext(nmod_poly_t delta, nmod_poly_t xi, slong n,
		const nmod_poly_t modulus, const mp_limb_t z) {

	if (n == 1) {
		nmod_poly_set(delta, delta_init_trivial);
		nmod_poly_set(xi, xi_init_trivial);

		return;
	}

	slong z_degree = 0;
	nmod_poly_t temp_delta;
	nmod_poly_t temp_xi;
	nmod_poly_init(temp_delta, modulus->mod.n);
	nmod_poly_init(temp_xi, modulus->mod.n);

	if (n % 2 == 0) {

		compute_semi_trace_trivial_ext(delta, xi, n / 2, modulus, z);
		nmod_poly_set(temp_delta, delta);
		nmod_poly_set(temp_xi, xi);
		z_degree = nmod_poly_degree(modulus) - n / 2;

	} else {

		compute_semi_trace_trivial_ext(delta, xi, n - 1, modulus, z);
		nmod_poly_set(temp_delta, delta_init_trivial);
		nmod_poly_set(temp_xi, xi_init_trivial);
		z_degree = nmod_poly_degree(modulus) - 1;

	}

	// compute delta
	nmod_poly_compose_mod(delta, delta, temp_xi, modulus);
	nmod_poly_scalar_mul_nmod(delta, delta, nmod_pow_ui(z, z_degree, modulus->mod));
	nmod_poly_add(delta, delta, temp_delta);

	//compute xi
	nmod_poly_compose_mod(xi, xi, temp_xi, modulus);

	nmod_poly_clear(temp_delta);
	nmod_poly_clear(temp_xi);
}

/**
 * Computes the value $\delta_r = a + z^{r - 1}\sigma(a) + z^{r - 2}\sigma^2(a) + \cdots + 
 * z^{1}\sigma^{r - 1}(a)$ where $a \in \mathbb{F}_p[x][z]$, and $r$ is the degree of
 * the extension {@code ctx}.
 */
void FFIsomorphism::compute_semi_trace_small_ext(fq_nmod_poly_t theta, const fq_nmod_poly_t a, const fq_nmod_ctx_t ctx,
		const fq_nmod_poly_t modulus) {

	fq_nmod_t xi;
	fq_nmod_init(xi, ctx);

	fq_nmod_poly_set(delta_init, a, ctx);
	compute_semi_trace_small_ext(theta, xi, fq_nmod_ctx_degree(ctx), ctx, modulus);

	fq_nmod_clear(xi, ctx);
}

void FFIsomorphism::compute_semi_trace_trivial_ext(nmod_poly_t theta, const nmod_poly_t a, const nmod_poly_t modulus,
		mp_limb_t z) {

	// use naive linear algebra for low-degree moduli
	if (nmod_poly_degree(modulus) < TRACE_THRESHOLD) {
		slong rows = nmod_poly_degree(modulus);
		nmod_mat_t frob_auto;
		nmod_mat_init(frob_auto, rows, rows, modulus->mod.n);

		nmod_poly_t temp;
		nmod_poly_init(temp, modulus->mod.n);
		nmod_poly_set(temp, xi_init_trivial);

		nmod_mat_entry(frob_auto, 0, 0) = 1;
		for (slong i = 1; i < rows; i++) {
			for (slong j = 0; j < rows; j++) {
				nmod_mat_entry(frob_auto, j, i) = nmod_poly_get_coeff_ui(temp, j);
			}

			nmod_poly_mulmod(temp, temp, xi_init_trivial, modulus);
		}

		for (slong i = 0; i < rows; i++)
			nmod_mat_entry(frob_auto, i, i) = nmod_sub(nmod_mat_entry(frob_auto, i, i), z, modulus->mod);

		nmod_mat_nullspace(frob_auto, frob_auto);

		nmod_poly_zero(theta);
		for (slong i = 0; i < rows; i++)
			nmod_poly_set_coeff_ui(theta, i, nmod_mat_entry(frob_auto, i, 0));

		nmod_mat_clear(frob_auto);
		nmod_poly_clear(temp);

		return;
	}

	nmod_poly_t xi;
	nmod_poly_init(xi, modulus->mod.n);
	nmod_poly_set(delta_init_trivial, a);

	compute_semi_trace_trivial_ext(theta, xi, nmod_poly_degree(modulus), modulus, z);

	nmod_poly_clear(xi);
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
 * Computes the value $\theta = \alpha + z^{r - 1}\alpha^{p} + z^{r - 2}\alpha^{p^2} + \cdots + 
 * z\alpha^{p^{(r - 1)}$ where $\alpha \in \mathbb{F}_p[x]$, and $r$ is the degree of
 * the extension {@code ctx}.
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
 * $\alpha^{p^i}$ for all $i = 0, \dots, r - 1$. The algorithm used is Algorithm 3.1
 * form Von Zur Gathen and Shoup, 1992. More precisely, this method implements the
 * case {q = p, R = ctx, t = p, m = r - 1} of the algorithm in the paper.
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
	fq_nmod_t factor;
	fq_nmod_init(factor, ctx_1);

	slong r = fq_nmod_ctx_degree(ctx_1);
	slong p = factor->mod.n;
	NModCyclotomicPoly nModCyclotomicPoly;
	nModCyclotomicPoly.single_irred_factor(factor, r, p);
	fq_nmod_ctx_init_modulus(cyclotomic_ctx, factor, "z");

	// build the modulus
	convert(modulus, cyclotomic_ctx->modulus, ctx_1);

	nmod_poly_clear(factor);
}

mp_limb_t FFIsomorphism::compute_cyclotomic_root(const slong r, const mp_limb_t p) {
	nmod_poly_t factor;
	nmod_poly_init(factor, p);

	NModCyclotomicPoly nModCyclotomicPoly;
	nModCyclotomicPoly.single_irred_factor(factor, r, p);
	mp_limb_t root = p - nmod_poly_get_coeff_ui(factor, 0);

	nmod_poly_clear(factor);

	return root;
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

	fq_nmod_t alpha;
	fq_nmod_init(alpha, ctx);

	// computing xi_init = x^p
	nmod_poly_set_coeff_ui(alpha, 1, 1);
	fq_nmod_pow_ui(xi_init, alpha, ctx->modulus->mod.n, ctx);

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

	// try alpha = x first
	nmod_poly_set_coeff_ui(alpha, 1, 1);
	compute_semi_trace_large_ext(theta, alpha, ctx, modulus);

	// if the semi trace of x is zero then we try random cases
	while (fq_nmod_poly_is_zero(theta, ctx)) {
		fq_nmod_randtest_not_zero(alpha, state, ctx);
		compute_semi_trace_large_ext(theta, alpha, ctx, modulus);
	}

	fq_nmod_clear(alpha, ctx);
	flint_randclear(state);
}

void FFIsomorphism::compute_semi_trace(nmod_poly_t theta, const nmod_poly_t modulus, const mp_limb_t z) {
	nmod_poly_t alpha;
	nmod_poly_init(alpha, modulus->mod.n);

	flint_rand_t state;
	flint_randinit(state);

	// set xi_init_trivial to x^p
	nmod_poly_set_coeff_ui(alpha, 1, 1);
	nmod_poly_powmod_ui_binexp(xi_init_trivial, alpha, modulus->mod.n, modulus);

	nmod_poly_zero(alpha);
	nmod_poly_zero(theta);

	while (nmod_poly_is_zero(theta)) {
		nmod_poly_randtest(alpha, state, nmod_poly_degree(modulus));
		compute_semi_trace_trivial_ext(theta, alpha, modulus, z);
	}
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

void FFIsomorphism::compute_extension_isomorphism(nmod_poly_t f, nmod_poly_t f_image) {
	nmod_poly_t temp;
	nmod_poly_init(temp, f->mod.n);

	slong r = fq_nmod_ctx_degree(ctx_1);
	mp_limb_t cyclo_root = compute_cyclotomic_root(r, f->mod.n);


	compute_semi_trace(f, ctx_1->modulus, cyclo_root);
	compute_semi_trace(f_image, ctx_2->modulus, cyclo_root);

	// compute the middle isomorphism
	nmod_poly_powmod_ui_binexp(temp, f, r, ctx_1->modulus);
	mp_limb_t eta1 = nmod_poly_get_coeff_ui(temp, 0);
	nmod_poly_powmod_ui_binexp(temp, f_image, r, ctx_2->modulus);
	mp_limb_t eta2 = nmod_poly_get_coeff_ui(temp, 0);
	mp_limb_t c = nmod_div(eta1, eta2, f->mod);

	CyclotomicExtRthRoot cyclotomicExtRthRoot;
	mp_limb_t root = cyclotomicExtRthRoot.compute_rth_root(c, r, f->mod.n);

	nmod_poly_scalar_mul_nmod(f_image, f_image, root);

	nmod_poly_clear(temp);
}

void FFIsomorphism::compute_generators(nmod_poly_t g1, nmod_poly_t g2) {

	// check for the trivial cyclotomic extension case
	Util util;
	slong degree = util.compute_multiplicative_order(g1->mod.n, fq_nmod_ctx_degree(ctx_1));
	if (degree == 1) {
		nmod_poly_t f;
		nmod_poly_t f_image;
		nmod_poly_init(f, g1->mod.n);
		nmod_poly_init(f_image, g1->mod.n);

		compute_extension_isomorphism(f, f_image);
		nmod_poly_set(g1, f);
		nmod_poly_set(g2, f_image);

		nmod_poly_clear(f);
		nmod_poly_clear(f_image);

		return;
	}

	fq_nmod_poly_t f;
	fq_nmod_poly_t f_image;

	fq_nmod_poly_init(f, ctx_1);
	fq_nmod_poly_init(f_image, ctx_2);

	compute_extension_isomorphism(f, f_image);

	fq_nmod_poly_get_coeff(g1, f, 0, ctx_1);
	fq_nmod_poly_get_coeff(g2, f_image, 0, ctx_2);

	fq_nmod_poly_clear(f, ctx_1);
	fq_nmod_poly_clear(f_image, ctx_2);
}

void FFIsomorphism::build_isomorphism(const nmod_poly_t g1, const nmod_poly_t g2) {
	fq_nmod_t x;
	fq_nmod_init(x, ctx_1);
	nmod_poly_set_coeff_ui(x, 1, 1);

	FFIsomBaseChange ffIsomBaseChange;
	ffIsomBaseChange.change_basis(x_image, g1, x, ctx_1->modulus);
	nmod_poly_compose_mod(x_image, x_image, g2, ctx_2->modulus);

	fq_nmod_clear(x, ctx_2);
}

void FFIsomorphism::compute_isom_matrix() {
	slong rows = nmod_poly_degree(ctx_2->modulus);
	nmod_poly_t temp;
	nmod_poly_init(temp, ctx_2->mod.n);
	nmod_poly_set(temp, x_image);

	nmod_mat_entry(isom_mat, 0, 0) = 1;

	for (slong i = 1; i < rows; i++) {
		for (slong j = 0; j < rows; j++) {
			nmod_mat_entry(isom_mat, j, i) = nmod_poly_get_coeff_ui(temp, j);
		}

		nmod_poly_mulmod(temp, temp, x_image, ctx_2->modulus);
	}

	nmod_poly_clear(temp);
}

void FFIsomorphism::get_x_image(nmod_poly_t x_image) {
	nmod_poly_set(x_image, this->x_image);
}

void FFIsomorphism::compute_image_using_modcomp(nmod_poly_t image, const nmod_poly_t f) {
	nmod_poly_compose_mod(image, f, x_image, ctx_2->modulus);
}

void FFIsomorphism::compute_image_using_matrix(nmod_poly_t image, const nmod_poly_t f) {
	nmod_mat_t f_vec;
	nmod_mat_t result_vec;
	slong rows = nmod_poly_degree(ctx_2->modulus);

	nmod_mat_init(f_vec, rows, 1, f->mod.n);
	nmod_mat_init(result_vec, rows, 1, f->mod.n);

	for (slong i = 0; i < rows; i++)
		nmod_mat_entry(f_vec, i, 0) = nmod_poly_get_coeff_ui(f, i);

	nmod_mat_mul(result_vec, isom_mat, f_vec);

	nmod_poly_zero(image);
	for (slong i = 0; i < rows; i++)
		nmod_poly_set_coeff_ui(image, i, nmod_mat_entry(result_vec, i, 0));

	nmod_mat_clear(f_vec);
	nmod_mat_clear(result_vec);
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
	fq_nmod_poly_init(delta_init, ctx_1);
	fq_nmod_init(xi_init, ctx_1);
	nmod_poly_init(delta_init_trivial, modulus1->mod.n);
	nmod_poly_init(xi_init_trivial, modulus1->mod.n);
	slong rows = nmod_poly_degree(ctx_2->modulus);
	nmod_mat_init(isom_mat, rows, rows, ctx_2->mod.n);

	nmod_poly_clear(tempf1);
	nmod_poly_clear(tempf2);
}

FFIsomorphism::~FFIsomorphism() {
	fq_nmod_poly_clear(delta_init, ctx_1);
	fq_nmod_clear(xi_init, ctx_1);
	nmod_poly_clear(delta_init_trivial);
	nmod_poly_clear(xi_init_trivial);
	fq_nmod_ctx_clear(ctx_1);
	fq_nmod_ctx_clear(ctx_2);
	nmod_poly_clear(x_image);
	nmod_mat_clear(isom_mat);
}
