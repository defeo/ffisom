/*
 * ff_isom_aux.h
 *
 *  Created on: Feb 11, 2014
 *      Author: javad
 */

#ifndef FF_ISOM_AUX_H_
#define FF_ISOM_AUX_H_

#include <flint/nmod_poly.h>

class FFIsomBaseChange {
    void monomial_to_dual(mp_limb_t *dual, const nmod_poly_t a, const nmod_poly_t modulus,
            const nmod_poly_t modulus_inv_rev);
    void dual_to_monomial(nmod_poly_t result, const mp_limb_t *dual, const nmod_poly_t modulus);

public:

    /**
     * Given {@code f} and {@code g} polynomials in $\mathbb{F}_p[X]/(modulus)$, computes
     * a polynomial $h \in \mathbb{F}_p[X]/(modulus)$ such that $h(f) = g$ if such polynomial exists. 
     * 
     * @param result	a polynomial h such that $h(f) = g$ if such polynomial exists
     */
    void change_basis(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t g,
            const nmod_poly_t modulus);
};

#endif /* FF_ISOM_AUX_H_ */
