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
	nmod_poly_t **R = bin_mat_init(result->mod.n);
	nmod_poly_t a;
	nmod_poly_t b;

	nmod_poly_init(a, result->mod.n);
	nmod_poly_init(b, result->mod.n);

	slong length = 2 * degree;

	for (slong i = 0; i < length; i++)
		nmod_poly_set_coeff_ui(a, i, sequence[i]);

	nmod_poly_set_coeff_ui(b, length, 1);

	halfgcd(R, b, a, degree);

	nmod_poly_reverse(result, R[1][1], nmod_poly_length(R[1][1]));
	nmod_poly_make_monic(result, result);

	bin_mat_clear(R);
	nmod_poly_clear(a);
	nmod_poly_clear(b);
}

/**
 * Computes the half gcd $R$ of $r_0, r_1$ where $r_0, r_1$. It is assumed that
 * $\deg(r_0) \ge \deg(r_1)$, and $0 \le k \le \deg(r_0)$.
 * 
 * @param R		The resulting half gcd
 * @return 		an integer h(k) depending on k such that the rows of R are the 
 * 				h-th and (h + 1)-th intermediate rows of the extended Euclidean algorithm 
 */
slong NmodMinPoly::halfgcd(nmod_poly_t **R, const nmod_poly_t r0, const nmod_poly_t r1, slong k) {

	slong temp = check_halfgcd_condition(R, r0, r1, k);
	if (temp >= 0)
		return temp;

	slong n0 = nmod_poly_degree(r0);
	slong n1 = nmod_poly_degree(r1);
	slong d = (slong) ceil((double) k / (double) 2);
	slong shift_index = 0;

	nmod_poly_t temp_r0;
	nmod_poly_t temp_r1;
	nmod_poly_t temp_r2;
	nmod_poly_t **tempR = bin_mat_init(r0->mod.n);

	nmod_poly_init(temp_r0, r0->mod.n);
	nmod_poly_init(temp_r1, r0->mod.n);
	nmod_poly_init(temp_r2, r0->mod.n);

	shift_index = n0 - 2 * d + 2;
	if (shift_index >= 0) {
		nmod_poly_shift_right(temp_r0, r0, shift_index);
		nmod_poly_shift_right(temp_r1, r1, shift_index);
	} else {
		nmod_poly_shift_left(temp_r0, r0, -shift_index);
		nmod_poly_shift_left(temp_r1, r1, -shift_index);
	}

	slong j = halfgcd(tempR, temp_r0, temp_r1, d - 1) + 1;
	slong delta = nmod_poly_degree(tempR[1][1]);

	shift_index = n0 - 2 * k;
	if (shift_index >= 0) {
		nmod_poly_shift_right(temp_r0, r0, shift_index);
		nmod_poly_shift_right(temp_r1, r1, shift_index);
	} else {
		nmod_poly_shift_left(temp_r0, r0, -shift_index);
		nmod_poly_shift_left(temp_r1, r1, -shift_index);
	}

	bin_mat_mul(tempR, temp_r0, temp_r1);
	n0 = nmod_poly_degree(temp_r0);
	n1 = nmod_poly_degree(temp_r1);

	if (nmod_poly_is_zero(temp_r1) || k < delta + n0 - n1) {
		bin_mat_set(R, tempR);

		nmod_poly_clear(temp_r0);
		nmod_poly_clear(temp_r1);
		nmod_poly_clear(temp_r2);
		bin_mat_clear(tempR);

		return j - 1;
	}

	nmod_poly_t q;
	nmod_poly_init(q, r0->mod.n);

	nmod_poly_divrem(q, temp_r2, temp_r0, temp_r1);

	nmod_poly_t **Q = bin_mat_init(r0->mod.n);
	nmod_poly_zero(Q[0][0]);
	nmod_poly_one(Q[0][1]);
	nmod_poly_one(Q[1][0]);
	nmod_poly_set(Q[1][1], q);
	nmod_poly_neg(Q[1][1], Q[1][1]);

	d = k - delta + n1 - n0;

	shift_index = n1 - 2 * d;
	if (shift_index >= 0) {
		nmod_poly_shift_right(temp_r1, temp_r1, shift_index);
		nmod_poly_shift_right(temp_r2, temp_r2, shift_index);
	} else {
		nmod_poly_shift_left(temp_r1, temp_r1, -shift_index);
		nmod_poly_shift_left(temp_r1, temp_r1, -shift_index);
	}

	nmod_poly_t **S = bin_mat_init(r0->mod.n);
	j = halfgcd(S, temp_r1, temp_r2, d) + j;

	bin_mat_mul(S, S, Q);
	bin_mat_mul(R, S, tempR);

	nmod_poly_clear(temp_r0);
	nmod_poly_clear(temp_r1);
	nmod_poly_clear(temp_r2);
	nmod_poly_clear(q);
	bin_mat_clear(tempR);
	bin_mat_clear(Q);
	bin_mat_clear(S);

	return j;
}

/**
 * Checks the base case conditions for the recursive halfgcd algorithm. 
 */
slong NmodMinPoly::check_halfgcd_condition(nmod_poly_t **R, const nmod_poly_t r0, const nmod_poly_t r1,
		slong k) {
	slong n0 = nmod_poly_degree(r0);
	slong n1 = nmod_poly_degree(r1);

	if (nmod_poly_is_zero(r1) || k < n0 - n1) {
		nmod_poly_one(R[0][0]);
		nmod_poly_zero(R[0][1]);
		nmod_poly_zero(R[1][0]);
		nmod_poly_one(R[1][1]);

		return 0;
	}

	if (k == 0 && n0 - n1 == 0) {
		nmod_poly_zero(R[0][0]);
		nmod_poly_one(R[0][1]);
		nmod_poly_one(R[1][0]);

		mp_limb_t temp0;
		mp_limb_t temp1;

		temp0 = nmod_poly_get_coeff_ui(r0, n0);
		temp1 = nmod_poly_get_coeff_ui(r1, n1);
		temp1 = nmod_inv(temp1, r1->mod);
		temp0 = nmod_mul(temp0, temp1, r1->mod);
		temp0 = nmod_neg(temp0, r1->mod);
		nmod_poly_zero(R[1][1]);
		nmod_poly_set_coeff_ui(R[1][1], 0, temp0);

		return 1;
	}

	return -1;
}

/**
 * Computes the matrix-matrix product $A = B \times C$.
 */
void NmodMinPoly::bin_mat_mul(nmod_poly_t **A, nmod_poly_t **B, nmod_poly_t **C) {
	nmod_poly_t **temp = bin_mat_init(B[0][0]->mod.n);
	bin_mat_set(temp, C);
	bin_mat_mul(B, temp[0][0], temp[1][0]);
	bin_mat_mul(B, temp[0][1], temp[1][1]);
	bin_mat_set(A, temp);

	bin_mat_clear(temp);
}

/**
 * Frees the memory allocated to the matrix {@code A}
 */
void NmodMinPoly::bin_mat_clear(nmod_poly_t **A) {
	for (slong i = 0; i < 2; i++) {
		for (slong j = 0; j < 2; j++)
			nmod_poly_clear(A[i][j]);
		delete[] A[i];
	}

	delete[] A;
}

/**
 * Allocates memory for a 2x2 matrix and initializes entries to zero. 
 */
nmod_poly_t** NmodMinPoly::bin_mat_init(const mp_limb_t p) {
	nmod_poly_t **bin_mat = new nmod_poly_t*[2];

	for (slong i = 0; i < 2; i++) {
		bin_mat[i] = new nmod_poly_t[2];

		for (slong j = 0; j < 2; j++)
			nmod_poly_init(bin_mat[i][j], p);
	}

	return bin_mat;
}

/**
 * Computes the matrix-vector product [f, g]^T = R * [f, g]^T.
 */
void NmodMinPoly::bin_mat_mul(nmod_poly_t **R, nmod_poly_t f, nmod_poly_t g) {
	nmod_poly_t temp1;
	nmod_poly_t temp2;
	nmod_poly_init(temp1, f->mod.n);
	nmod_poly_init(temp2, f->mod.n);

	nmod_poly_mul(temp1, R[0][0], f);
	nmod_poly_mul(temp2, R[0][1], g);
	nmod_poly_add(temp1, temp1, temp2);

	nmod_poly_mul(temp2, R[1][0], f);
	// don't need f anymore
	nmod_poly_set(f, temp1);
	nmod_poly_mul(temp1, R[1][1], g);
	nmod_poly_add(temp1, temp1, temp2);
	nmod_poly_set(g, temp1);

	nmod_poly_clear(temp1);
	nmod_poly_clear(temp2);
}

/**
 * Sets {@code A} to a compy of {@code B}
 */
void NmodMinPoly::bin_mat_set(nmod_poly_t **A, nmod_poly_t **B) {
	for (slong i = 0; i < 2; i++)
		for (slong j = 0; j < 2; j++)
			nmod_poly_set(A[i][j], B[i][j]);
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
