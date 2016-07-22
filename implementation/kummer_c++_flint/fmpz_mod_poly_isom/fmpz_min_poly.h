/*
 * fmpz_min_poly.h
 *
 *  Created on: Jan 22, 2014
 *      Author: javad
 */

#ifndef FMPZ_MIN_POLY_H_
#define FMPZ_MIN_POLY_H_

#include <flint/fmpz_vec.h>
#include <flint/fmpz_mod_poly.h>

class NmodMinPoly {

public:
	void inner_product(fmpz_t result, const fmpz_t *a, const fmpz_mod_poly_t f);
	void transposed_mul(fmpz_mod_poly_t result, const fmpz_mod_poly_t a, const fmpz_mod_poly_t b, slong n);
	void transposed_rem(fmpz_mod_poly_t result, const fmpz_mod_poly_t a, const fmpz_mod_poly_t r,
			const fmpz_mod_poly_t alpha,
			slong n);
	void transposed_mulmod(fmpz_t *result, const fmpz_t *r, const fmpz_mod_poly_t b, const fmpz_mod_poly_t a,
			const fmpz_mod_poly_t alpha);

	slong halfgcd(fmpz_mod_poly_t **r, const fmpz_mod_poly_t r0, const fmpz_mod_poly_t r1, slong k);
	slong check_halfgcd_condition(fmpz_mod_poly_t **R, const fmpz_mod_poly_t r0, const fmpz_mod_poly_t r1, slong k);
	fmpz_mod_poly_t **bin_mat_init(const fmpz_t p);
	void bin_mat_mul(fmpz_mod_poly_t **R, fmpz_mod_poly_t f, fmpz_mod_poly_t g);
	void bin_mat_set(fmpz_mod_poly_t **A, fmpz_mod_poly_t **B);
	void bin_mat_mul(fmpz_mod_poly_t **A, fmpz_mod_poly_t **B, fmpz_mod_poly_t **C);
	void bin_mat_clear(fmpz_mod_poly_t **A);

	void minimal_polynomial(fmpz_mod_poly_t result, const fmpz_t *sequence, slong length);


	void project_powers(fmpz_t *result, const fmpz_t *a, slong l, const fmpz_mod_poly_t h,
			const fmpz_mod_poly_t modulus, const fmpz_mod_poly_t modulus_inv_rev);
	void minimal_polynomial(fmpz_mod_poly_t result, const fmpz_mod_poly_t f, const fmpz_mod_poly_t modulus);
	void minimal_polynomial(fmpz_mod_poly_t result, const fmpz_mod_poly_t f, const fmpz_mod_poly_t modulus,
			const fmpz_mod_poly_t modulus_inv_rev);
};

#endif /* FMPZ_MIN_POLY_H_ */
