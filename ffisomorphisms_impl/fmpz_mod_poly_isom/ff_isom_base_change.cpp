/*
 * ff_isom_aux.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: javad
 */

#include "ff_isom_base_change.h"
#include "fmpz_min_poly.h"
#include <iostream>

using namespace std;


/**
 * Computes the dual representation of the given polynomial {@code a} modulo {@code modulus}
 * 
 * @param	modulus_inv_rev = 1 / rev(modulus, m + 1) mod x^m where m = deg(modulus)
 */
void FFIsomBaseChange::monomial_to_dual(fmpz_t *dual, const fmpz_mod_poly_t a, const fmpz_mod_poly_t modulus,
		const fmpz_mod_poly_t modulus_inv_rev) {
	fmpz_mod_poly_t temp;
	fmpz_mod_poly_init(temp, &modulus->p);

	slong m = fmpz_mod_poly_degree(modulus);

	fmpz_mod_poly_derivative(temp, modulus);
	fmpz_mod_poly_mulmod(temp, temp, a, modulus);
	fmpz_mod_poly_reverse(temp, temp, m);
	fmpz_mod_poly_mullow(temp, temp, modulus_inv_rev, m);

	for (slong i = 0; i < m; i++)
		fmpz_mod_poly_get_coeff_fmpz(dual[i], temp, i);

	fmpz_mod_poly_clear(temp);
}
/**
 * Computes the representation of the given dual {@code dual}
 * in monomial basis modulo {@code modulus} 
 */
void FFIsomBaseChange::dual_to_monomial(fmpz_mod_poly_t result, const fmpz_t *dual, const fmpz_mod_poly_t modulus) {
	fmpz_mod_poly_t temp1;
	fmpz_mod_poly_t temp2;
	fmpz_mod_poly_init(temp1, &modulus->p);
	fmpz_mod_poly_init(temp2, &modulus->p);

	slong m = fmpz_mod_poly_degree(modulus);
	for (slong i = 0; i < m; i++)
		fmpz_mod_poly_set_coeff_fmpz(temp1, i, dual[i]);
	fmpz_mod_poly_reverse(temp2, modulus, m + 1);
	fmpz_mod_poly_mullow(temp1, temp1, temp2, m);

	fmpz_mod_poly_reverse(temp2, temp1, m);

	// compute 1 / modulus'
	fmpz_mod_poly_derivative(temp1, modulus);
	fmpz_mod_poly_invmod(temp1, temp1, modulus);

	fmpz_mod_poly_mulmod(result, temp1, temp2, modulus);

	fmpz_mod_poly_clear(temp1);
	fmpz_mod_poly_clear(temp2);
}

void FFIsomBaseChange::change_basis(fmpz_mod_poly_t result, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g,
		const fmpz_mod_poly_t modulus) {
	fmpz_mod_poly_t min_poly;
	fmpz_mod_poly_init(min_poly, &f->p);

	fmpz_mod_poly_t modulus_inv_rev1;
	fmpz_mod_poly_t modulus_inv_rev2;
	fmpz_mod_poly_init(modulus_inv_rev1, &modulus->p);
	fmpz_mod_poly_init(modulus_inv_rev2, &modulus->p);

	slong degree = fmpz_mod_poly_degree(modulus);
	// compute 1 / rev(modulus, degree + 1) mod x^degree
	fmpz_mod_poly_reverse(modulus_inv_rev1, modulus, degree + 1);
	fmpz_mod_poly_inv_series_newton(modulus_inv_rev1, modulus_inv_rev1, degree);
	// compute 1 / rev(modulus, degree + 1) mod x^{degree - 1}
	fmpz_mod_poly_reverse(modulus_inv_rev2, modulus, degree + 1);
	fmpz_mod_poly_inv_series_newton(modulus_inv_rev2, modulus_inv_rev2, degree - 1);
	
	// compute the minimal polynomial of f
	NmodMinPoly fmpzMinPoly;
	fmpzMinPoly.minimal_polynomial(min_poly, f, modulus, modulus_inv_rev2);
	
	fmpz_t *dual = new fmpz_t[degree];
	for (slong i = 0; i < degree; i++)
		fmpz_init(dual[i]);

	// compute the dual basis image of x mod modulus
	monomial_to_dual(dual, g, modulus, modulus_inv_rev1);

	// compute the power projection <dual, f>
	fmpzMinPoly.project_powers(dual, dual, degree, f, modulus, modulus_inv_rev2);

	// compute result such that result(f) = x
	dual_to_monomial(result, dual, min_poly);

	fmpz_mod_poly_clear(min_poly);
	fmpz_mod_poly_clear(modulus_inv_rev1);
	fmpz_mod_poly_clear(modulus_inv_rev2);
	for (slong i = 0; i < degree; i++)
		fmpz_clear(dual[i]);
	delete[] dual;
}

