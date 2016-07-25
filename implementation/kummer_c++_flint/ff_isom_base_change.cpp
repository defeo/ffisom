/*
 * ff_isom_aux.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: javad
 */

#include "ff_isom_base_change.h"
#include "nmod_min_poly.h"
#include <iostream>

using namespace std;

/**
 * Computes the dual representation of the given polynomial {@code a} modulo {@code modulus}
 * 
 * @param	modulus_inv_rev = 1 / rev(modulus, m + 1) mod x^m where m = deg(modulus)
 */
void FFIsomBaseChange::monomial_to_dual(mp_limb_t *dual, const nmod_poly_t a, const nmod_poly_t modulus,
		const nmod_poly_t modulus_inv_rev) {
	nmod_poly_t temp;
	nmod_poly_init(temp, modulus->mod.n);

	slong m = nmod_poly_degree(modulus);

	nmod_poly_derivative(temp, modulus);
	nmod_poly_mulmod(temp, temp, a, modulus);
	nmod_poly_reverse(temp, temp, m);
	nmod_poly_mullow(temp, temp, modulus_inv_rev, m);

	for (slong i = 0; i < m; i++)
		dual[i] = nmod_poly_get_coeff_ui(temp, i);

	nmod_poly_clear(temp);
}

/**
 * Computes the representation of the given dual {@code dual}
 * in monomial basis modulo {@code modulus} 
 */
void FFIsomBaseChange::dual_to_monomial(nmod_poly_t result, const mp_limb_t *dual, const nmod_poly_t modulus) {
	nmod_poly_t temp1;
	nmod_poly_t temp2;
	nmod_poly_init(temp1, modulus->mod.n);
	nmod_poly_init(temp2, modulus->mod.n);

	slong m = nmod_poly_degree(modulus);
	for (slong i = 0; i < m; i++)
		nmod_poly_set_coeff_ui(temp1, i, dual[i]);

	nmod_poly_reverse(temp2, modulus, m + 1);
	nmod_poly_mullow(temp1, temp1, temp2, m);

	nmod_poly_reverse(temp2, temp1, m);

	// compute 1 / modulus'
	nmod_poly_derivative(temp1, modulus);
	nmod_poly_invmod(temp1, temp1, modulus);

	nmod_poly_mulmod(result, temp1, temp2, modulus);

	nmod_poly_clear(temp1);
	nmod_poly_clear(temp2);
}

void FFIsomBaseChange::change_basis(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t g,
		const nmod_poly_t modulus) {
	nmod_poly_t min_poly;
	nmod_poly_init(min_poly, f->mod.n);

	nmod_poly_t modulus_inv_rev1;
	nmod_poly_t modulus_inv_rev2;
	nmod_poly_init(modulus_inv_rev1, modulus->mod.n);
	nmod_poly_init(modulus_inv_rev2, modulus->mod.n);

	slong degree = nmod_poly_degree(modulus);
	// compute 1 / rev(modulus, degree + 1) mod x^degree
	nmod_poly_reverse(modulus_inv_rev1, modulus, degree + 1);
	nmod_poly_inv_series_newton(modulus_inv_rev1, modulus_inv_rev1, degree);
	// compute 1 / rev(modulus, degree + 1) mod x^{degree - 1}
	nmod_poly_reverse(modulus_inv_rev2, modulus, degree + 1);
	nmod_poly_inv_series_newton(modulus_inv_rev2, modulus_inv_rev2, degree - 1);

	// compute the minimal polynomial of f
	NmodMinPoly nmodMinPoly;
	nmodMinPoly.minimal_polynomial(min_poly, f, modulus, modulus_inv_rev2);

	mp_limb_t *dual = new mp_limb_t[degree];
	for (slong i = 0; i < degree; i++)
		dual[i] = 0;

	// compute the dual basis image of x mod modulus
	monomial_to_dual(dual, g, modulus, modulus_inv_rev1);

	// compute the power projection <dual, f>
	nmodMinPoly.project_powers(dual, dual, degree, f, modulus, modulus_inv_rev2);

	// compute result such that result(f) = x
	dual_to_monomial(result, dual, min_poly);

	nmod_poly_clear(min_poly);
	nmod_poly_clear(modulus_inv_rev1);
	nmod_poly_clear(modulus_inv_rev2);
	_nmod_vec_clear(dual);
}

