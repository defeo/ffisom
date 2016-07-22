/* 
 * File:   ff_isom_artin_schreier.h
 * Author: javad
 *
 * Created on December 19, 2015, 12:31 PM
 */

#ifndef FF_ISOM_ARTIN_SCHREIER_H
#define FF_ISOM_ARTIN_SCHREIER_H

#include <flint/nmod_poly.h>

class FFIsomArtinSchreier {
    nmod_poly_t modulus1;
    nmod_poly_t modulus2;
    
    void compute_hilbert_90_expression(nmod_poly_t xi, nmod_poly_t beta_a,
	    nmod_poly_t beta_theta, nmod_poly_t alpha, const nmod_poly_t xi_init,
	    const nmod_poly_t a, const nmod_poly_t theta, const nmod_poly_t modulus,
	    slong n);
    void compute_hilbert_90_solution(nmod_poly_t result, const nmod_poly_t a,
	    const nmod_poly_t xi_init, const nmod_poly_t modulus);
    void compute_generator(nmod_poly_t result, const nmod_poly_t modulus, slong exponent);

public:

    /**
     * @param f1 Defining modulus of the first extension
     * @param f2 Defining modulus of the second extension
     */
    FFIsomArtinSchreier(const nmod_poly_t f1, const nmod_poly_t f2);

    /**
     * Computes generators g1 of ctx_1, and g2 of ctx_2 such that
     * h: ctx_1 --> ctx_2
     *       g1 --> g2
     * is an isomorphism  
     * @param g1
     * @param g2
     */
    void compute_generators(nmod_poly_t g1, nmod_poly_t g2);

    ~FFIsomArtinSchreier();
};

#endif /* FF_ISOM_ARTIN_SCHREIER_H */

