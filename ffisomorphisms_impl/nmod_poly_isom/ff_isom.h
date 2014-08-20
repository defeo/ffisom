/*
 * ffIsomorphism.h
 *
 *  Created on: Dec 3, 2013
 *      Author: javad
 */

#ifndef FFISOMORPHISM_H_
#define FFISOMORPHISM_H_

#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>

class FFIsomorphism {

	fq_nmod_ctx_t ctx_1;
	fq_nmod_ctx_t ctx_2;

	fq_nmod_poly_t delta_init;
	fq_nmod_t xi_init;

	void compute_semi_trace_small_ext(fq_nmod_poly_t delta, fq_nmod_t xi, slong n, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus);
	void compute_semi_trace_small_ext(fq_nmod_poly_t theta, const fq_nmod_poly_t a, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus);
	void compute_semi_trace_large_ext(fq_nmod_poly_t theta, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus);
	void compute_semi_trace(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus);
	void compute_xi(fq_nmod_t xi, const fq_nmod_t old_xi, const fq_nmod_ctx_t ctx);
	void compute_delta(fq_nmod_poly_t delta, const fq_nmod_t xi, slong z_degree, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus);
	void compute_xi_init(const fq_nmod_ctx_t ctx, slong s);
	void iterated_frobenius(fq_nmod_t *result, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx, slong s);
	void compute_extension_isomorphism(fq_nmod_poly_t f, fq_nmod_poly_t f_image);
	void compute_middle_isomorphism(fq_nmod_t c, const fq_nmod_poly_t theta_a, const fq_nmod_poly_t theta_b, const fq_nmod_poly_t modulus,
			const fq_nmod_ctx_t cyclotomic_ctx);
	void build_cyclotomic_extension(fq_nmod_poly_t modulus, fq_nmod_ctx_t cyclotomic_ctx);
	void convert(fq_nmod_t result, const fq_nmod_poly_t value, const fq_nmod_ctx_t ctx);
	void convert(fq_nmod_poly_t result, const fq_nmod_t value, const fq_nmod_ctx_t ctx);
	void build_isomorphism();

public:

	FFIsomorphism(const nmod_poly_t f1, const nmod_poly_t f2);
	nmod_poly_t x_image;
	void compute_image(nmod_poly_t image, const nmod_poly_t f);
	~FFIsomorphism();
};

#endif /* FFISOMORPHISM_H_ */
