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
 * Computes the minimal polynomial of the sequence {@code sequence}.
 * 
 * @param result	the minimal polynomial of degree <= degree
 * @param sequence	a sequence of length >= {@code 2 * degree}
 */
void NmodMinPoly::minimal_polynomial(nmod_poly_t result, const mp_limb_t *sequence, slong degree) {
	slong length = 2 * degree;

	slong slength = length;
	while(slength > 0 && sequence[slength-1] == 0)
		slength--;

	mp_limb_t *xpow;
	xpow = (mp_limb_t *) malloc((length+1)*sizeof(mp_limb_t));
	for (slong i = 0; i < length; i++)
		xpow[i] = 0;
	xpow[length] = 1;

	mp_limb_t *A;
	A = (mp_limb_t *) malloc((length+1)*sizeof(mp_limb_t));
	slong lenA;
	mp_limb_t *B;
	B = (mp_limb_t *) malloc(slength*sizeof(mp_limb_t));
	slong lenB;

	mp_limb_t *Minv[4];
	for (slong i = 0; i < 4; i++)
		Minv[i] = (mp_limb_t *) malloc((length+1)*sizeof(mp_limb_t));
	slong lenMinv[4];

	_nmod_poly_hgcd(Minv, lenMinv, A, &lenA, B, &lenB, xpow, length+1, sequence, slength, result->mod);

	for (slong i = 0; i < lenMinv[0]; i++)
		nmod_poly_set_coeff_ui(result, lenMinv[0]-1-i, Minv[0][i]);
	nmod_poly_make_monic(result, result);

	free(xpow);
	free(A);
	free(B);
	for (slong i = 0; i < 4; i++)
		free(Minv[i]);
}

/*
 * Code for simple extension.
 */

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

/*
 * Same as above, different interface.
 */
void NmodMinPoly::minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const fq_nmod_ctx_t ctx) {

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
	mp_limb_t *sequence = _nmod_vec_init(2 * degree);

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
	mp_limb_t *temp_a = _nmod_vec_init(degree);
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
		nmod_poly_mulmod_preinv(h_powers[i], h_powers[i - 1], h, modulus, modulus_inv_rev);

	slong base = 0;
	for (slong i = 0; i < m; i++) {

		for (slong j = 0; (j < k) && (base + j < l); j++)
			result[base + j] = inner_product(temp_a, h_powers[j]);

		// todo: preconditioning on h_powers[k]
		transposed_mulmod(temp_a, temp_a, h_powers[k], modulus, modulus_inv_rev);
		base += k;
	}

	_nmod_vec_clear(temp_a);
	for (slong i = 0; i <= k; i++)
		nmod_poly_clear(h_powers[i]);
	delete[] h_powers;
}

/*
 * Same as below with different interface and preconditioning on a.
 * deg(a) < n
 * deg(b) < n
 * deg(ctx) = n
 */
void NmodMinPoly::transposed_mulmod_prerem(nmod_poly_t result, const nmod_poly_t arem, const nmod_poly_t b, const fq_nmod_ctx_t ctx) {

	slong n = fq_nmod_ctx_degree(ctx);

	transposed_mul(result, arem, b, n - 1);
}

/**
 * Computes the transposed modular product of {@code a} and {@code b} modulo {@code mod}: b°a.
 * 
 * @param a			a vector of length deg(mod)
 * @param b			a polynomial of degree < deg(mod)
 * @param mod			the modulus
 * @param mod_inv_rev		1/rev(mod, n) mod x^{n - 1} where n = deg(mod)
 */
void NmodMinPoly::transposed_mulmod(mp_limb_t *result, const mp_limb_t *a, const nmod_poly_t b, const nmod_poly_t mod,
		const nmod_poly_t mod_rev_inv) {

	slong n = nmod_poly_degree(mod);

	nmod_poly_t temp;
	nmod_poly_init2(temp, mod->mod.n, n);

	for (slong i = 0; i < n; i++)
		nmod_poly_set_coeff_ui(temp, i, a[i]);

	transposed_mulmod(temp, temp, b, mod, mod_rev_inv);

	for (slong i = 0; i < n; i++)
		result[i] = nmod_poly_get_coeff_ui(temp, i);

	nmod_poly_clear(temp);
}

/*
 * Same as above with different interface.
 * deg(a) < n
 * deg(b) < n
 * deg(ctx) = n
 */
void NmodMinPoly::transposed_mulmod(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t b, const fq_nmod_ctx_t ctx) {

	transposed_mulmod(result, a, b, ctx->modulus, ctx->inv);
}

/*
 * Same as above with different interface.
 * deg(a) < n
 * deg(b) < n
 * deg(ctx) = n
 */
void NmodMinPoly::transposed_mulmod(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t b, const nmod_poly_t mod, const nmod_poly_t mod_rev_inv) {

	slong m = nmod_poly_degree(b);
	if (m == -1) {
		nmod_poly_zero(result);
		return;
	}

	slong n = nmod_poly_degree(mod);

	nmod_poly_t temp;
	nmod_poly_init(temp, a->mod.n);

	transposed_rem(temp, a, mod, mod_rev_inv, m + n - 1);
	transposed_mul(result, temp, b, n - 1);

	nmod_poly_clear(temp);
}

/**
 * Computes the transposed product of {@code a} and {@code b}.
 * It is assumed that deg(a) <= deg(b) + m.
 * 
 * @param result	the transposed product of degree <= m
 */
void NmodMinPoly::transposed_mul(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t b,
		slong m) {

	slong n = nmod_poly_degree(b);
	if (n == -1) {
		nmod_poly_zero(result);
		return;
	}

	nmod_poly_t temp;
	nmod_poly_init(temp, b->mod.n);

	nmod_poly_reverse(temp, b, n + 1);
	nmod_poly_mullow(temp, a, temp, m + n + 1);
	nmod_poly_shift_right(temp, temp, n);
	nmod_poly_set(result, temp);

	nmod_poly_clear(temp);
}

/**
 * Computes the transposed remainder of {@code b} and {@code mod}.
 * 
 * @param a			the remainder of degree < n where n = deg(mod)
 * @param mod			the modulus of degree n
 * @param mod_inv_rev		1/rev(mod, n) mod x^{n - 1} where n = deg(mod)
 * @param result	the transposed remainder of degree <= m
 */
void NmodMinPoly::transposed_rem(nmod_poly_t result, const nmod_poly_t a,
		const nmod_poly_t mod,
		const nmod_poly_t mod_inv_rev,
		slong m) {

	nmod_poly_t temp;
	nmod_poly_init(temp, mod->mod.n);

	slong n = nmod_poly_degree(mod);

	transposed_mul(temp, a, mod, m - n);
	nmod_poly_mullow(temp, temp, mod_inv_rev, m - n + 1);
	nmod_poly_neg(temp, temp);
	nmod_poly_shift_left(temp, temp, n);
	nmod_poly_add(temp, a, temp);
	nmod_poly_set(result, temp);

	nmod_poly_clear(temp);
}

/**
 * Computes the inner product of {@code a} and {@code b} 
 */
mp_limb_t NmodMinPoly::inner_product(const nmod_poly_t a, const nmod_poly_t b) {
	slong len = FLINT_MIN(nmod_poly_degree(a), nmod_poly_degree(b)) + 1;
	slong nlimbs = _nmod_vec_dot_bound_limbs(len, b->mod);
	return _nmod_vec_dot(a->coeffs, b->coeffs, len, b->mod, nlimbs);
}

/**
 * Computes the inner product of {@code a} and the coefficient vector of {@code f} 
 */
mp_limb_t NmodMinPoly::inner_product(const mp_limb_t *a, const nmod_poly_t f) {
	slong len = nmod_poly_degree(f)+1;
	slong nlimbs = _nmod_vec_dot_bound_limbs(len, f->mod);
	return _nmod_vec_dot(a, f->coeffs, len, f->mod, nlimbs);
}


/*
 * Code for simple double extension.
 * Used for irreducible polynomial construction.
 */

void NmodMinPoly::minimal_polynomial(nmod_poly_t result, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx, slong m, const fq_nmod_t xi) {
	mp_limb_t *sequence = new mp_limb_t[2 * m];

	nmod_poly_t w;
	nmod_poly_init(w, ctx->modulus->mod.n);
	nmod_poly_one(w);

	fq_nmod_poly_t v;
	fq_nmod_poly_init(v, ctx);
	fq_nmod_poly_one(v, ctx);

	project_powers(sequence, w, v, 2*m, f, ctx, m, xi);
	minimal_polynomial(result, sequence, m);

	nmod_poly_clear(w);
	fq_nmod_poly_clear(v, ctx);

	_nmod_vec_clear(sequence);
}

/**
 * Power projection over a double extension of fields:
 * (F_p[x]/(irred_n))[y]/(y^m-xi)
 *
 * interface using polys to ease manipulation
 */
void NmodMinPoly::project_powers(mp_limb_t *result, const nmod_poly_t w, const fq_nmod_poly_t v, slong l, const fq_nmod_poly_t h,
		const fq_nmod_ctx_t ctx, slong m, const fq_nmod_t xi) {
	slong k = n_sqrt(l);
	slong lk = (slong) ceil((double) l / (double) k);

	// a temp for v to do computations without changing v
	fq_nmod_poly_t temp_v;
	fq_nmod_poly_init(temp_v, ctx);
	fq_nmod_poly_set(temp_v, v, ctx);
	// temp to store a coeff of v
	nmod_poly_t coeff_v;
	nmod_poly_init_preinv(coeff_v, w->mod.n, w->mod.ninv);

	// temp for transposed modular product
	nmod_poly_t *temp_w = new nmod_poly_t[m];
	for (slong i = 0; i < m; i++)
		nmod_poly_init_preinv(temp_w[i], w->mod.n, w->mod.ninv);

	// temp for modular product
	fq_nmod_poly_t high_part;
	fq_nmod_poly_init(high_part, ctx);

	// initials h_powers
	fq_nmod_poly_t *h_powers = new fq_nmod_poly_t[k + 1];
	for (slong i = 0; i <= k; i++)
		fq_nmod_poly_init(h_powers[i], ctx);
	// compute 1, h, h^2, ... h^k
	fq_nmod_poly_one(h_powers[0], ctx);
	fq_nmod_poly_set(h_powers[1], h, ctx);
	for (slong i = 2; i <= k; i++) {
		//multiply
		fq_nmod_poly_mul(h_powers[i], h_powers[i - 1], h, ctx);
		//reduce mod x^m - xi
		fq_nmod_poly_shift_right(high_part, h_powers[i], m, ctx);
		fq_nmod_poly_scalar_mul_fq_nmod(high_part, high_part, xi, ctx);
		fq_nmod_poly_truncate(h_powers[i], m, ctx);
		fq_nmod_poly_add(h_powers[i], h_powers[i], high_part, ctx);
	}
	// temp for coeff of h_powers
	nmod_poly_t coeff_h_power;
	nmod_poly_init_preinv(coeff_h_power, w->mod.n, w->mod.ninv);

	slong base = 0;
	for (slong i = 0; i < lk; i++) {
		for (slong t = 0; t < m; t++) {
			// accessing directly should temp_v->coeffs should work
			// even when the degree became too low
			fq_nmod_poly_get_coeff(coeff_v, temp_v, t, ctx);
			// todo: preconditioning on w
			transposed_mulmod(temp_w[t], w, coeff_v, ctx);
		}
	
		for (slong j = 0; (j < k) && (base + j < l); j++) {
			result[base + j] = 0;
			// we could avoid (fake) modular reduction in the sum
			for (slong t = 0; t < m; t++) {
				fq_nmod_poly_get_coeff(coeff_h_power, h_powers[j], t, ctx);
				result[base + j] = nmod_add(result[base + j], inner_product(temp_w[t], coeff_h_power), w->mod);
			}
		}

		// transposed modular multiplication in the double ext
		// todo: preconditioning on h_powers[k]
		transposed_mulmod(temp_v, temp_v, h_powers[k], ctx, m, xi);

		base += k;
	}

	fq_nmod_poly_clear(temp_v, ctx);
	nmod_poly_clear(coeff_v);

	for (slong i = 0; i < m; i++)
		nmod_poly_clear(temp_w[i]);
	delete [] temp_w;

	fq_nmod_poly_clear(high_part, ctx);
	for (slong i = 0; i <= k; i++)
		fq_nmod_poly_clear(h_powers[i], ctx);
	delete[] h_powers;
	nmod_poly_clear(coeff_h_power);
}

/*
 * Transposed modular multiplication for double extension.
 */
void NmodMinPoly::transposed_mulmod(fq_nmod_poly_t result, const fq_nmod_poly_t a, const fq_nmod_poly_t b, const fq_nmod_ctx_t ctx, slong m, const fq_nmod_t xi) {

	slong n = fq_nmod_poly_degree(b, ctx);
	if (n == -1) {
		fq_nmod_poly_zero(result, ctx);
		return;
	}

	fq_nmod_poly_t temp;
	fq_nmod_poly_init(temp, ctx);

	transposed_rem(temp, a, ctx, m, xi, n + m - 1);
	transposed_mul(result, temp, b, ctx, m - 1);

	fq_nmod_poly_clear(temp, ctx);
}

/**
 * Computes the transposed product of {@code a} and {@code b}.
 * It is assumed that deg(a) <= deg(b) + m.
 * 
 * @param result	the transposed product of degree <= m
 */
void NmodMinPoly::transposed_mul(fq_nmod_poly_t result, const fq_nmod_poly_t a, const fq_nmod_poly_t b, const fq_nmod_ctx_t ctx,
		slong m) {

	slong n = fq_nmod_poly_degree(b, ctx);
	if (n == -1) {
		fq_nmod_poly_zero(result, ctx);
		return;
	}

	fq_nmod_poly_t temp;
	fq_nmod_poly_init(temp, ctx);

	fq_nmod_poly_reverse(temp, b, n + 1, ctx);
	fq_nmod_poly_mullow(temp, a, temp, m + n + 1, ctx);
	fq_nmod_poly_shift_right(result, temp, n, ctx);

	fq_nmod_poly_clear(temp, ctx);
}

/**
 * Computes the transposed remainder of {@code b} and {@code x^m - xi}.
 * 
 * @param a			the remainder of degree < m
 * @param result	the transposed remainder of degree <= n
 */
void NmodMinPoly::transposed_rem(fq_nmod_poly_t result, const fq_nmod_poly_t a,
		const fq_nmod_ctx_t ctx, slong m, const fq_nmod_t xi,
		slong n) {
	// m and n switched

	fq_nmod_poly_t temp;
	fq_nmod_poly_init(temp, ctx);

        // transposed mul/reversed middle product with (x^m - xi)
	// only take -xi into account
	fq_nmod_poly_set(temp, a, ctx);
        fq_nmod_poly_truncate(temp, n - m + 2, ctx);
	fq_nmod_poly_scalar_mul_fq_nmod(temp, temp, xi, ctx);
	// the inverse of the reverse of (x^m - xi) is 1 + O(x^m)
	// no need for a mullow
	fq_nmod_poly_shift_left(temp, temp, m, ctx);
	fq_nmod_poly_add(result, a, temp, ctx);

	fq_nmod_poly_clear(temp, ctx);
}

