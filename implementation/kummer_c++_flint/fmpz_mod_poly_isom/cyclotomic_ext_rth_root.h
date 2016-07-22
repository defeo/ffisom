/*
 * fpz_rth_root.h
 *
 *  Created on: Dec 12, 2013
 *      Author: javad
 */

#ifndef FPZ_RTH_ROOT_H_
#define FPZ_RTH_ROOT_H_

#include <flint/fq_poly.h>
#include <flint/fq.h>


/**
 * This class is for computing an $r$-th root in the $r$-th cyclotomic extension of a finite
 * field $\mathbb{F}_p$ for arbitrary primes $r, p$. 
 */
class CyclotomicExtRthRoot {
	
	fq_poly_t beta_init;
	fq_poly_t xi_init;
	fq_t a;
	const fq_ctx_struct *ctx;
	slong r;
	
	void compute_trace(fq_poly_t beta, fq_poly_t xi, slong n, const fq_poly_t modulus);
	void compute_trace(fq_poly_t trace, const fq_poly_t alpha, const fq_poly_t modulus);
	void compute_xi(fq_poly_t xi, const fq_poly_t old_xi, slong zeta_degree);
	void compute_beta(fq_poly_t beta, const fq_poly_t xi, slong zeta_degree, const fq_poly_t modulus);
	void compute_beta_coeffs_small_ext(fq_poly_t beta, slong zeta_degree);
	void compute_beta_coeffs_large_ext(fq_poly_t beta, slong zeta_degree);
	void compute_beta_coeffs(fq_poly_t beta, slong zeta_degree);
	void compute_simple_compose(fq_t f, slong zeta_degree);
	void compute_simple_compose(fq_poly_t f, slong x_power, const fq_poly_t modulus);
	void compute_initials(const fq_poly_t alpha, const fq_poly_t modulus);
	void compute_rth_root_from_factor(fq_t root, const fq_poly_t factor);
	void compute_rth_root(fq_t root, fq_poly_t f);
	
public:
	
	/**
	 * Computes an $r$-th root of {@code a} in the cyclotomic field {@code ctx}.
	 * The field {@code ctx} is assumed to be a quotient $\mathbb{F}_p[Z] / (g(Z))$ 
	 * where $g(Z)$ is an irreducible factor of the $r$-th cyclotomic polynomial
	 * over $\mathbb{F}_p$. This method does not check for $r$-th powers, and assumes
	 * {@code a} an $r$-th power. 
	 * 
	 * @param root	the computed $r$-th root.
	 * @param a		an $r$-th power in the field {@code ctx}
	 * @param r		the order of root to be taken by
	 * @param ctx	a representation of the cyclotomic field
	 */
	void compute_rth_root(fq_t root, const fq_t a, slong r, const fq_ctx_t ctx);
	
};


#endif /* FPZ_RTH_ROOT_H_ */
