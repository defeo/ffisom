/*
 * ff_isom_aux.h
 *
 *  Created on: Feb 11, 2014
 *      Author: javad
 */

#ifndef FF_ISOM_AUX_H_
#define FF_ISOM_AUX_H_

#include <flint/fmpz_mod_poly.h>

class FFIsomBaseChange {

	void monomial_to_dual(fmpz_t *dual, const fmpz_mod_poly_t a, const fmpz_mod_poly_t modulus,
			const fmpz_mod_poly_t modulus_inv_rev);
	void dual_to_monomial(fmpz_mod_poly_t result, const fmpz_t *dual, const fmpz_mod_poly_t modulus);

public:

	/**
	 * Given {@code f} and {@code g} polynomials in $\mathbb{F}_p[X]/(modulus)$, computes
	 * a polynomial $h \in \mathbb{F}_p[X]/(modulus)$ such that $h(f) = g$ if such polynomial exists. 
	 * 
	 * @param result	the polynomial h such that $h(f) = g$ if such polynomial exists
	 */
	void change_basis(fmpz_mod_poly_t result, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g,
			const fmpz_mod_poly_t modulus);
};

#endif /* FF_ISOM_AUX_H_ */
