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

#include "nmod_poly_compose_mod.h"

enum {FORCE_LINALG, FORCE_MODCOMP, FORCE_COFACTOR, FORCE_ITERFROB, FORCE_MPE, FORCE_NONE};

class FFIsomPrimePower {
    slong ext_deg;
    mp_limb_t ext_char;
    fq_nmod_ctx_t ctx_1;
    fq_nmod_ctx_t ctx_2;

    slong cyclo_deg;
    nmod_poly_t cyclo_mod;
    fq_nmod_ctx_t cyclo_ctx;
    mp_limb_t cyclo_root;

    fq_nmod_poly_t delta_init;
    fq_nmod_t delta_init_trivial;
    fq_nmod_t xi_init;

    slong linalg_threshold;
    slong cofactor_threshold;
    slong iterfrob_threshold;
    slong mpe_threshold;

    void compute_semi_trace_trivial_linalg(fq_nmod_t theta, const fq_nmod_ctx_t ctx, mp_limb_t z);
    void compute_semi_trace_trivial_modcomp(fq_nmod_t theta, const fq_nmod_t a, const fq_nmod_ctx_t ctx, mp_limb_t z);
    void _compute_semi_trace_trivial_modcomp(fq_nmod_t delta, fq_nmod_t xi, slong n, const fq_nmod_ctx_t ctx, const mp_limb_t z);
    void compute_semi_trace_linalg_cyclo(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx);
    void lift_ht90_modexp_luca(fq_nmod_poly_t theta, const fq_nmod_t a, const fq_nmod_ctx_t ctx);
    void lift_ht90_modexp(fq_nmod_poly_t theta, const fq_nmod_t a0, const fq_nmod_ctx_t ctx);
    void lift_ht90_linalg(fq_nmod_poly_t theta, const fq_nmod_t a0, const nmod_mat_t frob_auto, const fq_nmod_ctx_t ctx);
     void compute_semi_trace_linalg(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx);
    void compute_semi_trace_cofactor_naive(fq_nmod_poly_t theta, const nmod_poly_t cofactor, const fq_nmod_ctx_t ctx);
    void compute_semi_trace_cofactor(fq_nmod_poly_t theta, const nmod_poly_t cofactor, const fq_nmod_ctx_t ctx);
    void compute_semi_trace_modcomp(fq_nmod_poly_t theta, const fq_nmod_poly_t a, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t cyclo_mod_lift);
    void _compute_semi_trace_modcomp(fq_nmod_poly_t delta, fq_nmod_t xi, slong n, 
				     const Nmod_poly_compose_mod & compose_xi_init, 
				     const fq_nmod_ctx_t ctx, const fq_nmod_poly_t cyclo_mod_lift);
    void compute_xi(fq_nmod_t xi, const fq_nmod_t old_xi, const fq_nmod_ctx_t ctx);
    void compute_delta(fq_nmod_poly_t delta, const fq_nmod_t xi, slong z_degree, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus);

    void compute_delta_and_xi(fq_nmod_poly_t delta, fq_nmod_t new_xi, const fq_nmod_t xi, slong z_degree, 
			      const Nmod_poly_compose_mod & compose_xi_init, 
			      const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus);

    void shift_delta(fq_nmod_poly_t delta, slong z_degree, const fq_nmod_ctx_t ctx);
    void compute_semi_trace_iterfrob_naive(fq_nmod_poly_t theta, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus);
    void compute_semi_trace_iterfrob(fq_nmod_poly_t theta, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t modulus);
    void iterated_frobenius(fq_nmod_t *result, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx, slong s);
    void compute_semi_trace(fq_nmod_t theta, const fq_nmod_ctx_t ctx, const mp_limb_t z);

    void compute_semi_trace(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t cyclo_mod_lift);
    void compute_semi_trace_all(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t cyclo_mod_lift);

    void compute_extension_isomorphism(fq_nmod_poly_t f, fq_nmod_poly_t f_image);
    void compute_middle_isomorphism(fq_nmod_t c, const fq_nmod_poly_t theta_a, const fq_nmod_poly_t theta_b, const fq_nmod_poly_t modulus);
    void convert(fq_nmod_t result, const fq_nmod_poly_t value, const fq_nmod_ctx_t ctx);
    void convert(fq_nmod_poly_t result, const fq_nmod_t value, const fq_nmod_ctx_t ctx);
    void compute_extension_isomorphism(nmod_poly_t f, nmod_poly_t f_image);
    void compute_cyclotomic_extension();
    void compute_cyclotomic_root();

public:

    /**
     * @param f1 Defining modulus of the first extension
     * @param f2 Defining modulus of the second extension
     */
    FFIsomPrimePower(const nmod_poly_t f1, const nmod_poly_t f2, slong force_algo = FORCE_NONE);

    /**
     * Computes generators g1 of ctx_1, and g2 of ctx_2 such that
     * h: ctx_1 --> ctx_2
     *       g1 --> g2
     * is an isomorphism  
     * @param g1
     * @param g2
     */
    void compute_generators_trivial(nmod_poly_t g1, nmod_poly_t g2);
    void compute_generators_nontriv(nmod_poly_t g1, nmod_poly_t g2);
    void compute_generators(nmod_poly_t g1, nmod_poly_t g2);

    ~FFIsomPrimePower();
};

#endif /* FFISOMORPHISM_H_ */
