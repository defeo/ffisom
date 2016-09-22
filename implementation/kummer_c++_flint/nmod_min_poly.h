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
#include <flint/fq_nmod_poly.h>

class NmodMinPoly {
public:
    mp_limb_t inner_product(const nmod_poly_t a, const nmod_poly_t b);
    mp_limb_t inner_product(const mp_limb_t *a, const nmod_poly_t f);
    void transposed_mul(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t b, slong n);
    void transposed_rem(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t r,
            const nmod_poly_t alpha,
            slong n);
    void transposed_mulmod_prerem(nmod_poly_t result, const nmod_poly_t arem, const nmod_poly_t b, const fq_nmod_ctx_t ctx);
    void transposed_mulmod(mp_limb_t *result, const mp_limb_t *r, const nmod_poly_t b, const nmod_poly_t a,
            const nmod_poly_t alpha);
    void transposed_mulmod(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t b, const fq_nmod_ctx_t ctx);
    void transposed_mulmod(nmod_poly_t result, const nmod_poly_t a, const nmod_poly_t b, const nmod_poly_t mod, const nmod_poly_t mod_rev_inv);
    void minimal_polynomial(nmod_poly_t result, const mp_limb_t *sequence, slong length);
    void project_powers(mp_limb_t *result, const mp_limb_t *a, slong l, const nmod_poly_t h,
            const nmod_poly_t modulus, const nmod_poly_t modulus_inv_rev);
    void minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const fq_nmod_ctx_t ctx);
    void minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t modulus);
    void minimal_polynomial(nmod_poly_t result, const nmod_poly_t f, const nmod_poly_t modulus,
            const nmod_poly_t modulus_inv_rev);

    void transposed_mul(fq_nmod_poly_t result, const fq_nmod_poly_t a, const fq_nmod_poly_t b, const fq_nmod_ctx_t ctx, slong m);
    void transposed_rem(fq_nmod_poly_t result, const fq_nmod_poly_t a, 
		const fq_nmod_ctx_t ctx, slong m, const fq_nmod_t xi,
		slong n);
    void transposed_mulmod(fq_nmod_poly_t result, const fq_nmod_poly_t a, const fq_nmod_poly_t b, const fq_nmod_ctx_t ctx, slong m, const fq_nmod_t xi);
    void project_powers(mp_limb_t *result, const nmod_poly_t w, const fq_nmod_poly_t v, slong l, const fq_nmod_poly_t h,
		const fq_nmod_ctx_t ctx, slong m, const fq_nmod_t xi);
    void minimal_polynomial(nmod_poly_t result, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx, slong m, const fq_nmod_t xi);
};

#endif /* NMOD_MIN_POLY_H_ */
