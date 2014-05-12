/*
 * ffIsomorphism.h
 *
 *  Created on: Dec 3, 2013
 *      Author: javad
 */

#ifndef FFISOMORPHISM_H_
#define FFISOMORPHISM_H_

#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>

class FFIsomorphism {

	fq_ctx_t ctx_1;
	fq_ctx_t ctx_2;
	fmpz_mod_poly_t x_image;

	fq_poly_t delta_init;
	fq_t xi_init;

	void compute_semi_trace_small_ext(fq_poly_t delta, fq_t xi, slong n, const fq_ctx_t ctx, const fq_poly_t modulus);
	void compute_semi_trace_small_ext(fq_poly_t theta, const fq_poly_t a, const fq_ctx_t ctx, const fq_poly_t modulus);
	void compute_semi_trace_large_ext(fq_poly_t theta, const fq_t alpha, const fq_ctx_t ctx, const fq_poly_t modulus);
	void compute_semi_trace(fq_poly_t theta, const fq_ctx_t ctx, const fq_poly_t modulus);
	void compute_xi(fq_t xi, const fq_t old_xi, const fq_ctx_t ctx);
	void compute_delta(fq_poly_t delta, const fq_t xi, slong z_degree, const fq_ctx_t ctx, const fq_poly_t modulus);
	void compute_xi_init(const fq_ctx_t ctx, slong s);
	void iterated_frobenius(fq_t *result, const fq_t alpha, const fq_ctx_t ctx, slong s);
	void compute_extension_isomorphism(fq_poly_t f, fq_poly_t f_image);
	void compute_middle_isomorphism(fq_t c, const fq_poly_t theta_a, const fq_poly_t theta_b, const fq_poly_t modulus,
			const fq_ctx_t cyclotomic_ctx);
	void build_cyclotomic_extension(fq_poly_t modulus, fq_ctx_t cyclotomic_ctx);
	void convert(fq_t result, const fq_poly_t value, const fq_ctx_t ctx);
	void convert(fq_poly_t result, const fq_t value, const fq_ctx_t ctx);
	void build_isomorphism();

public:

	FFIsomorphism(const fmpz_mod_poly_t f1, const fmpz_mod_poly_t f2);
	void compute_image(fmpz_mod_poly_t image, const fmpz_mod_poly_t f);
	~FFIsomorphism();
};

#endif /* FFISOMORPHISM_H_ */
