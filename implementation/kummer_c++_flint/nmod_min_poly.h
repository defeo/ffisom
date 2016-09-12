/*
 * fmpz_min_poly.h
 *
 *  Created on: Jan 22, 2014
 *      Author: javad
 */

#ifndef NMOD_MIN_POLY_H_
#define NMOD_MIN_POLY_H_

#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

class NmodMinPoly {
public:
    mp_limb_t inner_product(const mp_limb_t *a, const nmod_poly_t f);
    void transposed_mul(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t b, slong n);
    void transposed_rem(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t r,
            const nmod_poly_t alpha,
            slong n);
    void transposed_mulmod(mp_limb_t *result, const mp_limb_t *r, const nmod_poly_t b, const nmod_poly_t a,
            const nmod_poly_t alpha);
    slong halfgcd(nmod_poly_t **r, const nmod_poly_t r0, const nmod_poly_t r1, slong k);
    slong check_halfgcd_condition(nmod_poly_t **R, const nmod_poly_t r0, const nmod_poly_t r1, slong k);
    nmod_poly_t **bin_mat_init(const mp_limb_t p);
    void bin_mat_mul(nmod_poly_t **R, nmod_poly_t f, nmod_poly_t g);
    void bin_mat_set(nmod_poly_t **A, nmod_poly_t **B);
    void bin_mat_mul(nmod_poly_t **A, nmod_poly_t **B, nmod_poly_t **C);
    void bin_mat_clear(nmod_poly_t **A);
    void minimal_polynomial(nmod_poly_t result, const mp_limb_t *sequence, slong length);
    void project_powers(mp_limb_t *result, const mp_limb_t *a, slong l, const nmod_poly_t h,
            const nmod_poly_t modulus, const nmod_poly_t modulus_inv_rev);
    void minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const fq_nmod_t ctx);
    void minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t modulus);
    void minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t modulus,
            const nmod_poly_t modulus_inv_rev);
};

#endif /* NMOD_MIN_POLY_H_ */
