/*
 * fmpz_min_poly.cpp
 *
 *  Created on: Jan 22, 2014
 *      Author: javad
 */

#include "fmpz_min_poly.h"
#include <math.h>
#include <iostream>
#include <flint/fmpz_poly_mat.h>

using namespace std;

/**
 * Computes the minimal polynomial of {@code f} modulo {@code modulus}.
 * 
 * @param result	the minimal polynomial of {@code f}
 */

void NmodMinPoly::minimal_polynomial(fmpz_mod_poly_t result, const fmpz_mod_poly_t f, const fmpz_mod_poly_t modulus) {

	fmpz_mod_poly_t alpha;
	fmpz_mod_poly_init(alpha, &modulus->p);
	slong degree = fmpz_mod_poly_degree(modulus);
	// compute alpha = 1 / rev(degree, modulus) mod x^{degree - 1}
	fmpz_mod_poly_reverse(alpha, modulus, degree + 1);
	fmpz_mod_poly_inv_series_newton(alpha, alpha, degree - 1);

	minimal_polynomial(result, f, modulus, alpha);
	
	fmpz_mod_poly_clear(alpha);
}

/**
 * Computes the minimal polynomial of {@code f} modulo {@code modulus}.
 * 
 * @param result			the minimal polynomial of {@code f}
 * @param modulus_inv_rev	1 / rev(d, modulus) mod x^{d - 1} where d = deg(modulus)
 */
void NmodMinPoly::minimal_polynomial(fmpz_mod_poly_t result, const fmpz_mod_poly_t f, const fmpz_mod_poly_t modulus,
		const fmpz_mod_poly_t modulus_inv_rev) {
	fmpz_mod_poly_t g;
	fmpz_mod_poly_t temp_g;
	fmpz_mod_poly_t tau;

	fmpz_mod_poly_init(g, &f->p);
	fmpz_mod_poly_init(temp_g, &f->p);
	fmpz_mod_poly_init(tau, &f->p);

	fmpz_mod_poly_set_ui(g, 1);
	fmpz_mod_poly_set_ui(tau, 1);

	slong degree = fmpz_mod_poly_degree(modulus);
	fmpz_t *sequence = new fmpz_t[2 * degree];
	for (slong i = 0; i < 2 * degree; i++)
		fmpz_init(sequence[i]);

	flint_rand_t state;
	flint_randinit(state);

	slong l = 0;
	while (true) {
		
		for (slong i = 0; i < 2 * degree; i++)
			fmpz_randm(sequence[i], state, &f->p);

		transposed_mulmod(sequence, sequence, tau, modulus, modulus_inv_rev);
		l = degree - fmpz_mod_poly_degree(g);
		project_powers(sequence, sequence, l * 2, f, modulus, modulus_inv_rev);
		minimal_polynomial(temp_g, sequence, l);

		fmpz_mod_poly_mul(g, g, temp_g);
		if (fmpz_mod_poly_degree(g) == degree)
			break;

		fmpz_mod_poly_compose_mod(temp_g, temp_g, f, modulus);
		fmpz_mod_poly_mulmod(tau, tau, temp_g, modulus);
		if (fmpz_mod_poly_is_zero(tau))
			break;
	}

	fmpz_mod_poly_set(result, g);

	fmpz_mod_poly_clear(g);
	fmpz_mod_poly_clear(temp_g);
	fmpz_mod_poly_clear(tau);
	for (slong i = 0; i < 2 * degree; i++)
		fmpz_clear(sequence[i]);
	delete[] sequence;
	flint_randclear(state);
}

/**
 * Computes the minimal polynomial of the sequence {@code sequence}.
 * 
 * @param result	the minimal polynomial of degree <= degree
 * @param sequence	a sequence of length >= {@code 2 * degree}
 */
void NmodMinPoly::minimal_polynomial(fmpz_mod_poly_t result, const fmpz_t *sequence, slong degree) {
	fmpz_mod_poly_t **R = bin_mat_init(&result->p);
	fmpz_mod_poly_t a;
	fmpz_mod_poly_t b;

	fmpz_mod_poly_init(a, &result->p);
	fmpz_mod_poly_init(b, &result->p);
	
	slong length = 2 * degree;

	for (slong i = 0; i < length; i++)
		fmpz_mod_poly_set_coeff_fmpz(a, i, sequence[i]);
	
	fmpz_mod_poly_set_coeff_ui(b, length, 1);

	halfgcd(R, b, a, degree);

	fmpz_mod_poly_reverse(result, R[1][1], fmpz_mod_poly_length(R[1][1]));
	fmpz_mod_poly_make_monic(result, result);

	bin_mat_clear(R);
	fmpz_mod_poly_clear(a);
	fmpz_mod_poly_clear(b);
}

/**
 * Computes the half gcd $R$ of $r_0, r_1$ where $r_0, r_1$. It is assumed that
 * $\deg(r_0) \ge \deg(r_1)$, and $0 \le k \le \deg(r_0)$.
 * 
 * @param R		The resulting half gcd
 * @return 		an integer h(k) depending on k such that the rows of R are the 
 * 				h-th and (h + 1)-th intermediate rows of the extended Euclidean algorithm 
 */
slong NmodMinPoly::halfgcd(fmpz_mod_poly_t **R, const fmpz_mod_poly_t r0, const fmpz_mod_poly_t r1, slong k) {

	slong temp = check_halfgcd_condition(R, r0, r1, k);
	if (temp >= 0)
		return temp;

	slong n0 = fmpz_mod_poly_degree(r0);
	slong n1 = fmpz_mod_poly_degree(r1);
	slong d = (slong) ceil((double) k / (double) 2);
	slong shift_index = 0;

	fmpz_mod_poly_t temp_r0;
	fmpz_mod_poly_t temp_r1;
	fmpz_mod_poly_t temp_r2;
	fmpz_mod_poly_t **tempR = bin_mat_init(&r0->p);

	fmpz_mod_poly_init(temp_r0, &r0->p);
	fmpz_mod_poly_init(temp_r1, &r0->p);
	fmpz_mod_poly_init(temp_r2, &r0->p);

	shift_index = n0 - 2 * d + 2;
	if (shift_index >= 0) {
		fmpz_mod_poly_shift_right(temp_r0, r0, shift_index);
		fmpz_mod_poly_shift_right(temp_r1, r1, shift_index);
	} else {
		fmpz_mod_poly_shift_left(temp_r0, r0, -shift_index);
		fmpz_mod_poly_shift_left(temp_r1, r1, -shift_index);
	}

	slong j = halfgcd(tempR, temp_r0, temp_r1, d - 1) + 1;
	slong delta = fmpz_mod_poly_degree(tempR[1][1]);

	shift_index = n0 - 2 * k;
	if (shift_index >= 0) {
		fmpz_mod_poly_shift_right(temp_r0, r0, shift_index);
		fmpz_mod_poly_shift_right(temp_r1, r1, shift_index);
	} else {
		fmpz_mod_poly_shift_left(temp_r0, r0, -shift_index);
		fmpz_mod_poly_shift_left(temp_r1, r1, -shift_index);
	}

	bin_mat_mul(tempR, temp_r0, temp_r1);
	n0 = fmpz_mod_poly_degree(temp_r0);
	n1 = fmpz_mod_poly_degree(temp_r1);

	if (fmpz_mod_poly_is_zero(temp_r1) || k < delta + n0 - n1) {
		bin_mat_set(R, tempR);

		fmpz_mod_poly_clear(temp_r0);
		fmpz_mod_poly_clear(temp_r1);
		fmpz_mod_poly_clear(temp_r2);
		bin_mat_clear(tempR);

		return j - 1;
	}

	fmpz_mod_poly_t q;
	fmpz_mod_poly_init(q, &r0->p);

	fmpz_mod_poly_divrem(q, temp_r2, temp_r0, temp_r1);

	fmpz_mod_poly_t **Q = bin_mat_init(&r0->p);
	fmpz_mod_poly_set_ui(Q[0][0], 0);
	fmpz_mod_poly_set_ui(Q[0][1], 1);
	fmpz_mod_poly_set_ui(Q[1][0], 1);
	fmpz_mod_poly_set(Q[1][1], q);
	fmpz_mod_poly_neg(Q[1][1], Q[1][1]);

	d = k - delta + n1 - n0;

	shift_index = n1 - 2 * d;
	if (shift_index >= 0) {
		fmpz_mod_poly_shift_right(temp_r1, temp_r1, shift_index);
		fmpz_mod_poly_shift_right(temp_r2, temp_r2, shift_index);
	} else {
		fmpz_mod_poly_shift_left(temp_r1, temp_r1, -shift_index);
		fmpz_mod_poly_shift_left(temp_r1, temp_r1, -shift_index);
	}

	fmpz_mod_poly_t **S = bin_mat_init(&r0->p);
	j = halfgcd(S, temp_r1, temp_r2, d) + j;

	bin_mat_mul(S, S, Q);
	bin_mat_mul(R, S, tempR);

	fmpz_mod_poly_clear(temp_r0);
	fmpz_mod_poly_clear(temp_r1);
	fmpz_mod_poly_clear(temp_r2);
	fmpz_mod_poly_clear(q);
	bin_mat_clear(tempR);
	bin_mat_clear(Q);
	bin_mat_clear(S);

	return j;
}

/**
 * Checks the base case conditions for the recursive halfgcd algorithm. 
 */
slong NmodMinPoly::check_halfgcd_condition(fmpz_mod_poly_t **R, const fmpz_mod_poly_t r0, const fmpz_mod_poly_t r1,
slong k) {
	slong n0 = fmpz_mod_poly_degree(r0);
	slong n1 = fmpz_mod_poly_degree(r1);

	if (fmpz_mod_poly_is_zero(r1) || k < n0 - n1) {
		fmpz_mod_poly_set_ui(R[0][0], 1);
		fmpz_mod_poly_set_ui(R[0][1], 0);
		fmpz_mod_poly_set_ui(R[1][0], 0);
		fmpz_mod_poly_set_ui(R[1][1], 1);

		return 0;
	}

	if (k == 0 && n0 - n1 == 0) {
		fmpz_mod_poly_set_ui(R[0][0], 0);
		fmpz_mod_poly_set_ui(R[0][1], 1);
		fmpz_mod_poly_set_ui(R[1][0], 1);

		fmpz_t temp0;
		fmpz_t temp1;
		fmpz_init(temp0);
		fmpz_init(temp1);

		fmpz_mod_poly_get_coeff_fmpz(temp0, r0, n0);
		fmpz_mod_poly_get_coeff_fmpz(temp1, r1, n1);
		fmpz_invmod(temp1, temp1, &r1->p);
		fmpz_mul(temp0, temp0, temp1);
		fmpz_neg(temp0, temp0);
		fmpz_mod_poly_set_fmpz(R[1][1], temp0);

		fmpz_clear(temp0);
		fmpz_clear(temp1);

		return 1;
	}

	return -1;
}

/**
 * Computes the matrix-matrix product $A = B \times C$.
 */
void NmodMinPoly::bin_mat_mul(fmpz_mod_poly_t **A, fmpz_mod_poly_t **B, fmpz_mod_poly_t **C) {
	fmpz_mod_poly_t **temp = bin_mat_init(&B[0][0]->p);
	bin_mat_set(temp, C);
	bin_mat_mul(B, temp[0][0], temp[1][0]);
	bin_mat_mul(B, temp[0][1], temp[1][1]);
	bin_mat_set(A, temp);

	bin_mat_clear(temp);
}

/**
 * Frees the memory allocated to the matrix {@code A}
 */
void NmodMinPoly::bin_mat_clear(fmpz_mod_poly_t **A) {
	for (slong i = 0; i < 2; i++) {
		for (slong j = 0; j < 2; j++)
			fmpz_mod_poly_clear(A[i][j]);
		delete[] A[i];
	}

	delete[] A;
}

/**
 * Allocates memory for a 2x2 matrix and initializes entries to zero. 
 */
fmpz_mod_poly_t** NmodMinPoly::bin_mat_init(const fmpz_t p) {
	fmpz_mod_poly_t **bin_mat = new fmpz_mod_poly_t*[2];

	for (slong i = 0; i < 2; i++) {
		bin_mat[i] = new fmpz_mod_poly_t[2];

		for (slong j = 0; j < 2; j++)
			fmpz_mod_poly_init(bin_mat[i][j], p);
	}

	return bin_mat;
}

/**
 * Computes the matrix-vector product [f, g]^T = R * [f, g]^T.
 */
void NmodMinPoly::bin_mat_mul(fmpz_mod_poly_t **R, fmpz_mod_poly_t f, fmpz_mod_poly_t g) {
	fmpz_mod_poly_t temp1;
	fmpz_mod_poly_t temp2;
	fmpz_mod_poly_init(temp1, &f->p);
	fmpz_mod_poly_init(temp2, &f->p);

	fmpz_mod_poly_mul(temp1, R[0][0], f);
	fmpz_mod_poly_mul(temp2, R[0][1], g);
	fmpz_mod_poly_add(temp1, temp1, temp2);

	fmpz_mod_poly_mul(temp2, R[1][0], f);
	// don't need f anymore
	fmpz_mod_poly_set(f, temp1);
	fmpz_mod_poly_mul(temp1, R[1][1], g);
	fmpz_mod_poly_add(temp1, temp1, temp2);
	fmpz_mod_poly_set(g, temp1);

	fmpz_mod_poly_clear(temp1);
	fmpz_mod_poly_clear(temp2);
}

/**
 * Sets {@code A} to a compy of {@code B}
 */
void NmodMinPoly::bin_mat_set(fmpz_mod_poly_t **A, fmpz_mod_poly_t **B) {
	for (slong i = 0; i < 2; i++)
		for (slong j = 0; j < 2; j++)
			fmpz_mod_poly_set(A[i][j], B[i][j]);
}

/**
 * Computes the power projection:
 * \[\langle v, 1 \rangle, \langle v, h \rangle, \dots, \langle v, h^{l - 1} \rangle\]
 * where $\langle ., . \rangle$ is the inner product, considering $h^i$ as the vector
 * of its coefficients.
 * 
 * @param modulus_inv_rev	1 / rev(m + 1, modulus) mod x^{m - 1} where m = deg(modulus)
 */
void NmodMinPoly::project_powers(fmpz_t *result, const fmpz_t *a, slong l, const fmpz_mod_poly_t h,
		const fmpz_mod_poly_t modulus, const fmpz_mod_poly_t modulus_inv_rev) {
	slong k = n_sqrt(l);
	slong m = (slong) ceil((double) l / (double) k);
	slong degree = fmpz_mod_poly_degree(modulus);

	// a temp for a to do computations without changing a
	fmpz_t *temp_a = new fmpz_t[degree];
	for (slong i = 0; i < degree; i++)
		fmpz_init_set(temp_a[i], a[i]);

	fmpz_mod_poly_t *h_powers = new fmpz_mod_poly_t[k + 1];
	// initials h_powers
	for (slong i = 0; i <= k; i++)
		fmpz_mod_poly_init(h_powers[i], &h->p);

	// compute 1, h, h^2, ... h^k
	fmpz_mod_poly_set_ui(h_powers[0], 1);
	fmpz_mod_poly_set(h_powers[1], h);
	for (slong i = 2; i <= k; i++)
		fmpz_mod_poly_mulmod(h_powers[i], h_powers[i - 1], h, modulus);

	slong base = 0;
	for (slong i = 0; i < m; i++) {

		for (slong j = 0; (j < k) && (base + j < l); j++)
			inner_product(result[base + j], temp_a, h_powers[j]);

		transposed_mulmod(temp_a, temp_a, h_powers[k], modulus, modulus_inv_rev);
		base += k;
	}

	for (slong i = 0; i < degree; i++)
		fmpz_clear(temp_a[i]);
	delete[] temp_a;

	for (slong i = 0; i <= k; i++)
		fmpz_mod_poly_clear(h_powers[i]);
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
void NmodMinPoly::transposed_mulmod(fmpz_t *result, const fmpz_t *r, const fmpz_mod_poly_t b, const fmpz_mod_poly_t a,
		const fmpz_mod_poly_t alpha) {

	slong m = fmpz_mod_poly_degree(a);
	slong n = fmpz_mod_poly_degree(b);

	fmpz_mod_poly_t temp;
	fmpz_mod_poly_init2(temp, &a->p, m);

	for (slong i = 0; i < m; i++)
		fmpz_mod_poly_set_coeff_fmpz(temp, i, r[i]);

	transposed_rem(temp, temp, a, alpha, n + m - 1);
	transposed_mul(temp, temp, b, m - 1);

	for (slong i = 0; i < m; i++)
		fmpz_mod_poly_get_coeff_fmpz(result[i], temp, i);

	fmpz_mod_poly_clear(temp);
}

/**
 * Computes the transposed product of {@code c} and {@code a}. It is assumed that deg(a) = m
 * and deg(c) <= m + n.
 * 
 * @param result	the transposed product of degree <= n
 */
void NmodMinPoly::transposed_mul(fmpz_mod_poly_t result, const fmpz_mod_poly_t c, const fmpz_mod_poly_t a,
slong n) {

	slong m = fmpz_mod_poly_degree(a);
	fmpz_mod_poly_t temp;
	fmpz_mod_poly_init(temp, &a->p);

	fmpz_mod_poly_reverse(temp, a, m + 1);
	fmpz_mod_poly_mullow(temp, c, temp, m + n + 1);
	fmpz_mod_poly_shift_right(temp, temp, m);
	fmpz_mod_poly_set(result, temp);

	fmpz_mod_poly_clear(temp);
}
/**
 * Computes the transposed remainder of {@code r} and {@code a}.
 * 
 * @param r			the remainder of degree < m where m = deg(a)
 * @param a			the modulus of degree m
 * @param alpha		1/rev(a, m) mod x^{n - m + 1} where m = deg(a)
 * @param result	the transposed remainder
 */
void NmodMinPoly::transposed_rem(fmpz_mod_poly_t result, const fmpz_mod_poly_t r, const fmpz_mod_poly_t a,
		const fmpz_mod_poly_t alpha,
		slong n) {

	fmpz_mod_poly_t temp;
	fmpz_mod_poly_init(temp, &a->p);

	slong m = fmpz_mod_poly_degree(a);

	transposed_mul(temp, r, a, n - m);
	fmpz_mod_poly_mullow(temp, temp, alpha, n - m + 1);
	fmpz_mod_poly_shift_left(temp, temp, m);
	fmpz_mod_poly_neg(temp, temp);
	fmpz_mod_poly_add(temp, temp, r);
	fmpz_mod_poly_set(result, temp);

	fmpz_mod_poly_clear(temp);
}

/**
 * Computes the inner product of {@code a} and the coefficient vector of {@code f} 
 */
void NmodMinPoly::inner_product(fmpz_t result, const fmpz_t *a, const fmpz_mod_poly_t f) {
	fmpz_t temp;
	fmpz_init(temp);
	fmpz_zero(result);

	slong degree = fmpz_mod_poly_degree(f);

	for (slong i = 0; i <= degree; i++) {
		fmpz_mod_poly_get_coeff_fmpz(temp, f, i);
		fmpz_addmul(result, temp, a[i]);
	}
	fmpz_mod(result, result, &f->p);

	fmpz_clear(temp);
}

/////////////////////////////////// tests /////////////////////////////////////////////


//void test_power_projection() {
//	slong l = 34;
//	slong degree = 37;
//
//	fmpz_t p;
//	fmpz_init_set_ui(p, 9001);
//
//	fmpz_t *r = new fmpz_t[degree];
//	for (slong i = 0; i < degree; i++)
//		fmpz_init(r[i]);
//
//	fmpz_t *result = new fmpz_t[l];
//	for (slong i = 0; i < l; i++)
//		fmpz_init(result[i]);
//
//	flint_rand_t state;
//	flint_randinit(state);
//
//	for (slong i = 0; i < degree; i++)
//		fmpz_randm(r[i], state, p);
//
//	fmpz_mod_poly_t h, modulus, modulus_inv_rev;
//	fmpz_mod_poly_init(h, p);
//	fmpz_mod_poly_init(modulus, p);
//	fmpz_mod_poly_init(modulus_inv_rev, p);
//	
//	fmpz_mod_poly_randtest(h, state, degree);
//	fmpz_mod_poly_randtest(modulus, state, degree + 1);
//	fmpz_mod_poly_set_coeff_ui(modulus, degree, n_randtest(state));
//
//	fmpz_mod_poly_reverse(modulus_inv_rev, modulus, degree + 1);
//	fmpz_mod_poly_inv_series_newton(modulus_inv_rev, modulus_inv_rev, degree - 1);
//
//	FmpzMinPoly fmpzMinPoly;
//	fmpzMinPoly.project_powers(result, r, l, h, modulus, modulus_inv_rev);
//
//// naive ---------------------------------------------
//	fmpz_mod_poly_t *h_powers = new fmpz_mod_poly_t[l];
//	// initials h_powers
//	for (slong i = 0; i < l; i++)
//		fmpz_mod_poly_init(h_powers[i], &h->p);
//
//	// compute 1, h, h^2, ... h^(l - 1)
//	fmpz_mod_poly_set_ui(h_powers[0], 1);
//	fmpz_mod_poly_set(h_powers[1], h);
//	for (slong i = 2; i < l; i++)
//		fmpz_mod_poly_mulmod(h_powers[i], h_powers[i - 1], h, modulus);
//
//	fmpz_t *naive_result = new fmpz_t[l];
//	for (slong i = 0; i < l; i++)
//		fmpz_init(naive_result[i]);
//
//	for (slong i = 0; i < l; i++)
//		fmpzMinPoly.inner_product(naive_result[i], r, h_powers[i]);
//
//	_fmpz_vec_print(*naive_result, l);
//	cout << "\n";
//	_fmpz_vec_print(*result, l);
//	cout << "\n";
//
//	if (_fmpz_vec_equal(*result, *naive_result, l))
//		cout << "ok\n";
//	else
//		cout << "ooops\n";
//}
//
//
//
//void test_halfgcd() {
//	slong k = 10;
//	slong degree = 20;
//
//	fmpz_t p;
//	fmpz_init_set_ui(p, 9001);
//
//	flint_rand_t state;
//	flint_randinit(state);
//
//	fmpz_mod_poly_t f, g;
//	fmpz_mod_poly_init(f, p);
//	fmpz_mod_poly_init(g, p);
//	fmpz_mod_poly_randtest(f, state, degree);
//	fmpz_mod_poly_randtest(g, state, degree - 5);
//
//	FmpzMinPoly fmpzMinPoly;
//	fmpz_mod_poly_t **R = fmpzMinPoly.bin_mat_init(p);
//	k = fmpz_mod_poly_degree(f);
//	fmpzMinPoly.halfgcd(R, f, g, k);
//
//	fmpz_mod_poly_mul(R[1][0], R[0][0], f);
//	fmpz_mod_poly_mul(R[1][1], R[0][1], g);
//	fmpz_mod_poly_add(R[1][0], R[1][0], R[1][1]);
//	fmpz_t a;
//	fmpz_init(a);
//	fmpz_mod_poly_get_coeff_fmpz(a, R[1][0], fmpz_mod_poly_degree(R[1][0]));
//	fmpz_invmod(a, a, p);
//	fmpz_mod_poly_scalar_mul_fmpz(R[1][0], R[1][0], a);
//	fmpz_mod_poly_scalar_mul_fmpz(R[0][0], R[0][0], a);
//	fmpz_mod_poly_scalar_mul_fmpz(R[0][1], R[0][1], a);
//
//	fmpz_mod_poly_print_pretty(R[0][0], "x");
//	cout << "\n";
//	fmpz_mod_poly_print_pretty(R[0][1], "x");
//	cout << "\n";
//	fmpz_mod_poly_print_pretty(R[1][0], "x");
//	cout << "\n";
//
//	fmpz_mod_poly_t gcd, s, t;
//	fmpz_mod_poly_init(gcd, p);
//	fmpz_mod_poly_init(s, p);
//	fmpz_mod_poly_init(t, p);
//
//	fmpz_mod_poly_xgcd(gcd, s, t, f, g);
//	fmpz_mod_poly_print_pretty(s, "x");
//	cout << "\n";
//	fmpz_mod_poly_print_pretty(t, "x");
//	cout << "\n";
//}
