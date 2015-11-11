/* 
 * File:   nmod_poly_cyclotomy.h
 * Author: javad
 *
 * Created on November 8, 2015, 3:10 PM
 */

#ifndef NMOD_POLY_CYCLOTOMY_H
#define NMOD_POLY_CYCLOTOMY_H

#include <flint/nmod_poly.h>

class NModCyclotomicPoly {
    slong n;
    slong s;

    void compose(nmod_poly_t result, const nmod_poly_t f, slong n);
    void equal_degree_fact(nmod_poly_factor_t factors, const nmod_poly_t f, flint_rand_t state);
    void single_irred_factor(nmod_poly_t factor, const nmod_poly_t f, flint_rand_t state);
    void compute_trace(nmod_poly_t result, const nmod_poly_t g, slong i, const nmod_poly_t modulus);
    void compute_power(nmod_poly_t result, const nmod_poly_t g, slong i);
    void split(nmod_poly_t f1, nmod_poly_t f2, const nmod_poly_t f, flint_rand_t state);

public:

    /**
     * Constructs the p-th cyclotomic polynomial where {@code p} is prime. 
     * @param result
     * @param p
     */
    void construct_cyclo_prime_degree(nmod_poly_t result, slong p);

    /**
     * Constructs the p^j-th cyclotomic polynomial where {@code p} is prime. 
     * @param result
     * @param p
     * @param j
     */
    void construct_cyclo_prime_power_degree(nmod_poly_t result, slong p, slong j);

    /**
     * Constructs the n-th cyclotomic polynomial where {@code n} is a positive integer.
     * The algorithm is quasi-linear with complexity O(M(n)log(n)). See "von zur Gathen, 
     * Gerhard, Modern Computer Algebra, Chapter 14".
     * @param result
     * @param n
     */
    void construct_cyclo(nmod_poly_t result, slong n);

    /**
     * Obtains all irreducible factors of the n-th cyclotomic polynomial with 
     * coefficients modulo {@code modulus}. The implementation is based on a fast
     * variant of Cantor-Zassenhaus Equal-Degree-Factorization. See "Shoup, Fast 
     * construction of irreducible polynomials over finite fields, Section 4"
     * @param factors
     * @param n
     * @param moulus
     */
    void all_irred_factors(nmod_poly_factor_t factors, slong n, slong moulus);

    /**
     * Obtains a single irreducible factor of the n-th cyclotomic polynomial with 
     * coefficients modulo {@code modulus}. See {@link #all_irred_factors() all_irred_factors}
     * method.
     * @param factor
     * @param n
     * @param modulus
     */
    void single_irred_factor(nmod_poly_t factor, slong n, slong modulus);
};

#endif /* NMOD_POLY_CYCLOTOMY_H */

