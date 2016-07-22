/*
 * fpz_rth_root.h
 *
 *  Created on: Dec 12, 2013
 *      Author: javad
 */

#ifndef FPZ_RTH_ROOT_H_
#define FPZ_RTH_ROOT_H_

#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod.h>

/**
 * This class is for computing an $r$-th root in the $r$-th cyclotomic extension of a finite
 * field $\mathbb{F}_p$ for arbitrary primes $r, p$. 
 */
class CyclotomicExtRthRoot {
    fq_nmod_poly_t beta_init;
    fq_nmod_poly_t xi_init;
    fq_nmod_t a;
    const fq_nmod_ctx_struct *ctx;
    slong r;

    void compute_trace(fq_nmod_poly_t beta, fq_nmod_poly_t xi, slong n, const fq_nmod_poly_t modulus);
    void compute_trace(fq_nmod_poly_t trace, const fq_nmod_poly_t alpha, const fq_nmod_poly_t modulus);
    void compute_xi(fq_nmod_poly_t xi, const fq_nmod_poly_t old_xi, slong zeta_degree);
    void compute_beta(fq_nmod_poly_t beta, const fq_nmod_poly_t xi, slong zeta_degree, const fq_nmod_poly_t modulus);
    void compute_beta_coeffs_small_ext(fq_nmod_poly_t beta, slong zeta_degree);
    void compute_beta_coeffs_large_ext(fq_nmod_poly_t beta, slong zeta_degree);
    void compute_beta_coeffs(fq_nmod_poly_t beta, slong zeta_degree);
    void compute_simple_compose(fq_nmod_t f, slong zeta_degree);
    void compute_simple_compose(fq_nmod_poly_t f, slong x_power, const fq_nmod_poly_t modulus);
    void compute_initials(const fq_nmod_poly_t alpha, const fq_nmod_poly_t modulus);
    void compute_rth_root_from_factor(fq_nmod_t root, const fq_nmod_poly_t factor);
    void compute_rth_root(fq_nmod_t root, const fq_nmod_poly_t f);
    mp_limb_t compute_rth_root_from_factor(const mp_limb_t a, const nmod_poly_t factor);
    mp_limb_t compute_rth_root(const nmod_poly_t f);

public:

    /**
     * Computes an $r$-th root of {@code a} in the cyclotomic field {@code ctx}.
     * The field {@code ctx} is assumed to be a quotient $\mathbb{F}_p[Z] / (g(Z))$ 
     * where $g(Z)$ is an irreducible factor of the $r$-th cyclotomic polynomial
     * over $\mathbb{F}_p$. This method does not check for $r$-th powers, and assumes
     * {@code a} an $r$-th power. 
     * 
     * @param root  the computed $r$-th root.
     * @param a	    an $r$-th power in the field {@code ctx}
     * @param r	    the order of root to be taken by
     * @param ctx   a representation of the cyclotomic field
     */
    void compute_rth_root(fq_nmod_t root, const fq_nmod_t a, slong r, const fq_nmod_ctx_t ctx);
    mp_limb_t compute_rth_root(const mp_limb_t c, slong r, slong p);

};


#endif /* FPZ_RTH_ROOT_H_ */
