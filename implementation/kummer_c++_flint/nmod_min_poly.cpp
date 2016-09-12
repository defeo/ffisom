/*
 * fmpz_min_poly.cpp
 *
 *  Created on: Jan 22, 2014
 *      Author: javad
 */

#include "nmod_min_poly.h"
#include <math.h>
#include <iostream>
#include <flint/nmod_poly_mat.h>

using namespace std;

/**
 * Computes the minimal polynomial of {@code f} modulo {@code modulus}.
 * 
 * @param result	the minimal polynomial of {@code f}
 */

void NmodMinPoly::minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t modulus) {

	nmod_poly_t alpha;
	nmod_poly_init(alpha, modulus->mod.n);
	slong degree = nmod_poly_degree(modulus);
	// compute alpha = 1 / rev(degree, modulus) mod x^{degree - 1}
	nmod_poly_reverse(alpha, modulus, degree + 1);
	nmod_poly_inv_series_newton(alpha, alpha, degree - 1);

	minimal_polynomial(result, f, modulus, alpha);

	nmod_poly_clear(alpha);
}

void NmodMinPoly::minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const fq_nmod_t ctx) {

	minimal_polynomial(result, f, ctx->modulus, ctx->inv);
}

/**
 * Computes the minimal polynomial of {@code f} modulo {@code modulus}.
 * 
 * @param result			the minimal polynomial of {@code f}
 * @param modulus_inv_rev	1 / rev(d, modulus) mod x^{d - 1} where d = deg(modulus)
 */
void NmodMinPoly::minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t modulus,
		const nmod_poly_t modulus_inv_rev) {
	nmod_poly_t g;
	nmod_poly_t temp_g;
	nmod_poly_t tau;

	nmod_poly_init(g, f->mod.n);
	nmod_poly_init(temp_g, f->mod.n);
	nmod_poly_init(tau, f->mod.n);

	nmod_poly_one(g);
	nmod_poly_one(tau);

	slong degree = nmod_poly_degree(modulus);
	mp_limb_t *sequence = new mp_limb_t[2 * degree];

	flint_rand_t state;
	flint_randinit(state);

	slong l = 0;
	while (true) {

		_nmod_vec_randtest(sequence, state, 2 * degree, f->mod);

		transposed_mulmod(sequence, sequence, tau, modulus, modulus_inv_rev);
		l = degree - nmod_poly_degree(g);
		project_powers(sequence, sequence, l * 2, f, modulus, modulus_inv_rev);
		minimal_polynomial(temp_g, sequence, l);

		nmod_poly_mul(g, g, temp_g);
		if (nmod_poly_degree(g) == degree)
			break;

		nmod_poly_compose_mod(temp_g, temp_g, f, modulus);
		nmod_poly_mulmod(tau, tau, temp_g, modulus);
		if (nmod_poly_is_zero(tau))
			break;
	}

	nmod_poly_set(result, g);

	nmod_poly_clear(g);
	nmod_poly_clear(temp_g);
	nmod_poly_clear(tau);
	_nmod_vec_clear(sequence);
	flint_randclear(state);
}

/**
 * Computes the minimal polynomial of the sequence {@code sequence}.
 * 
 * @param result	the minimal polynomial of degree <= degree
 * @param sequence	a sequence of length >= {@code 2 * degree}
 */
void NmodMinPoly::minimal_polynomial(nmod_poly_t result, const mp_limb_t *sequence, slong degree) {
	slong length = 2 * degree;

	mp_limb_t *xpow;
	xpow = (mp_limb_t *) malloc((length+1)*sizeof(mp_limb_t));
	for (slong i = 0; i < length; i++)
		xpow[i] = 0;
	xpow[length] = 1;

	mp_limb_t *A;
	A = (mp_limb_t *) malloc((length+1)*sizeof(mp_limb_t));
	slong lenA;
	mp_limb_t *B;
	B = (mp_limb_t *) malloc(length*sizeof(mp_limb_t));
	slong lenB;

	mp_limb_t *Minv[4];
	for (slong i = 0; i < 4; i++)
		Minv[i] = (mp_limb_t *) malloc((length+1)*sizeof(mp_limb_t));
	slong lenMinv[4];

	_nmod_poly_hgcd(Minv, lenMinv, A, &lenA, B, &lenB, xpow, length+1, sequence, length, result->mod);

	for (slong i = 0; i < lenMinv[0]; i++)
		nmod_poly_set_coeff_ui(result, lenMinv[0]-1-i, Minv[0][i]);
	nmod_poly_make_monic(result, result);

	free(xpow);
	free(A);
	free(B);
	for (slong i = 0; i < 4; i++)
		free(Minv[i]);
}

/**
 * Computes the power projection:
 * \[\langle v, 1 \rangle, \langle v, h \rangle, \dots, \langle v, h^{l - 1} \rangle\]
 * where $\langle ., . \rangle$ is the inner product, considering $h^i$ as the vector
 * of its coefficients.
 * 
 * @param modulus_inv_rev	1 / rev(m + 1, modulus) mod x^{m - 1} where m = deg(modulus)
 */
void NmodMinPoly::project_powers(mp_limb_t *result, const mp_limb_t *a, slong l, const nmod_poly_t h,
		const nmod_poly_t modulus, const nmod_poly_t modulus_inv_rev) {
	slong k = n_sqrt(l);
	slong m = (slong) ceil((double) l / (double) k);
	slong degree = nmod_poly_degree(modulus);

	// a temp for a to do computations without changing a
	mp_limb_t *temp_a = new mp_limb_t[degree];
	for (slong i = 0; i < degree; i++)
		temp_a[i] = a[i];

	nmod_poly_t *h_powers = new nmod_poly_t[k + 1];
	// initials h_powers
	for (slong i = 0; i <= k; i++)
		nmod_poly_init(h_powers[i], h->mod.n);

	// compute 1, h, h^2, ... h^k
	nmod_poly_set_coeff_ui(h_powers[0], 0, 1);
	nmod_poly_set(h_powers[1], h);
	for (slong i = 2; i <= k; i++)
		nmod_poly_mulmod(h_powers[i], h_powers[i - 1], h, modulus);

	slong base = 0;
	for (slong i = 0; i < m; i++) {

		for (slong j = 0; (j < k) && (base + j < l); j++)
			result[base + j] = inner_product(temp_a, h_powers[j]);

		transposed_mulmod(temp_a, temp_a, h_powers[k], modulus, modulus_inv_rev);
		base += k;
	}

	_nmod_vec_clear(temp_a);
	for (slong i = 0; i <= k; i++)
		nmod_poly_clear(h_powers[i]);
	delete[] h_powers;
}

/**
 * Computes the transposed modular product of {@code r} and {@code b} modulo {@code a}.
 * 
 * @param r			a vector of length deg(a)
 * @param b			a polynomial of degree < deg(a)
 * @param a			the modulus
 * @param alpha		1/rev(a, m) mod x^{m - 1} where m = deg(a)
 */
void NmodMinPoly::transposed_mulmod(mp_limb_t *result, const mp_limb_t *r, const nmod_poly_t b, const nmod_poly_t a,
		const nmod_poly_t alpha) {

	slong m = nmod_poly_degree(a);
	slong n = nmod_poly_degree(b);

	nmod_poly_t temp;
	nmod_poly_init2(temp, a->mod.n, m);

	for (slong i = 0; i < m; i++)
		nmod_poly_set_coeff_ui(temp, i, r[i]);

	transposed_rem(temp, temp, a, alpha, n + m - 1);
	transposed_mul(temp, temp, b, m - 1);

	for (slong i = 0; i < m; i++)
		result[i] = nmod_poly_get_coeff_ui(temp, i);

	nmod_poly_clear(temp);
}

/**
 * Computes the transposed product of {@code c} and {@code a}. It is assumed that deg(a) = m
 * and deg(c) <= m + n.
 * 
 * @param result	the transposed product of degree <= n
 */
void NmodMinPoly::transposed_mul(nmod_poly_t result, const nmod_poly_t c, const nmod_poly_t a,
		slong n) {

	slong m = nmod_poly_degree(a);
	nmod_poly_t temp;
	nmod_poly_init(temp, a->mod.n);

	nmod_poly_reverse(temp, a, m + 1);
	nmod_poly_mullow(temp, c, temp, m + n + 1);
	nmod_poly_shift_right(temp, temp, m);
	nmod_poly_set(result, temp);

	nmod_poly_clear(temp);
}

/**
 * Computes the transposed remainder of {@code r} and {@code a}.
 * 
 * @param r			the remainder of degree < m where m = deg(a)
 * @param a			the modulus of degree m
 * @param alpha		1/rev(a, m) mod x^{n - m + 1} where m = deg(a)
 * @param result	the transposed remainder
 */
void NmodMinPoly::transposed_rem(nmod_poly_t result, const nmod_poly_t r, const nmod_poly_t a,
		const nmod_poly_t alpha,
		slong n) {

	nmod_poly_t temp;
	nmod_poly_init(temp, a->mod.n);

	slong m = nmod_poly_degree(a);

	transposed_mul(temp, r, a, n - m);
	nmod_poly_mullow(temp, temp, alpha, n - m + 1);
	nmod_poly_shift_left(temp, temp, m);
	nmod_poly_neg(temp, temp);
	nmod_poly_add(temp, temp, r);
	nmod_poly_set(result, temp);

	nmod_poly_clear(temp);
}

/**
 * Computes the inner product of {@code a} and the coefficient vector of {@code f} 
 */
mp_limb_t NmodMinPoly::inner_product(const mp_limb_t *a, const nmod_poly_t f) {
	mp_limb_t temp = 0;
	mp_limb_t result = 0;
	slong degree = nmod_poly_degree(f);

	for (slong i = 0; i <= degree; i++) {
		temp = nmod_poly_get_coeff_ui(f, i);
		temp = nmod_mul(temp, a[i], f->mod);
		result = nmod_add(result, temp, f->mod);
	}

	return result;
}
