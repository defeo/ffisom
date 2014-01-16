#ifndef NMOD_POLY_EXTRA_H
#define NMOD_POLY_EXTRA_H

#include <flint/nmod_poly.h>

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*              A few extra functionalities for Fp[x]                     */
/* 1. random polynomials                                                  */
/* 2. power series inverses needed for rem and trem                       */
/* 3. ratrecon and minimal polynomials                                    */
/* 4. transposed algorithms                                               */
/* 5. to and from Newtom sums                                             */
/* 6. misc. calculations mod M(x)                                         */
/* 7. geometric interpolation                                             */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/



/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* 1. random polynomials                                                  */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* random dense polynomial with at most len terms (so degree < len)       */
/*------------------------------------------------------------------------*/
void nmod_poly_rand_dense(nmod_poly_t poly, flint_rand_t state, long len);

/*------------------------------------------------------------------------*/
/* random dense monic polynomial of degree len                            */
/*------------------------------------------------------------------------*/
void nmod_poly_rand_dense_monic(nmod_poly_t poly, flint_rand_t state, long len);




/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* 2. power series inverses needed for rem and trem                       */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* res = 1 / reverse(P) mod x^k                                           */
/*------------------------------------------------------------------------*/
void nmod_poly_inverse_reverse(nmod_poly_t res, const nmod_poly_t P, const long k);

/*------------------------------------------------------------------------*/
/* res = 1 / reverse(P) mod x^{deg(P)-1}                                  */
/* (used for modular multiplication and its transpose)                    */
/*------------------------------------------------------------------------*/
void nmod_poly_inverse_reverse_main(nmod_poly_t res, const nmod_poly_t P);




/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* 3. ratrecon and minimal polynomials                                    */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* rational reconstruction                                                */
/* s must have degree <= 2d                                               */
/* output: deg(num) <= d, deg(den) <= d                                   */
/* BUGBUGBUGBUGBUG: sigma                                                 */
/*------------------------------------------------------------------------*/
void nmod_poly_ratrecon(nmod_poly_t num, nmod_poly_t den, const nmod_poly_t s, long d);

/*------------------------------------------------------------------------*/
/* minimal polynomial of degree at most d of the sequence val             */
/* val must have 2d+1 (or more) entries                                   */
/* (in reality, 2d would be enough but the code doesn't support it)       */
/*------------------------------------------------------------------------*/
void nmod_poly_minimal_polynomial_sequence(nmod_poly_t res, mp_srcptr val, long d);




/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* 4. transposed algorithms                                               */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* transposed product of ell by poly2                                     */
/* input: len(ell) = k+m, deg(poly2) <= m                                 */
/* output: len(res) = k                                                   */
/*------------------------------------------------------------------------*/
void nmod_poly_tmul(mp_ptr res, 
		    mp_srcptr ell, const nmod_poly_t poly2, const long m, const long k);

/*------------------------------------------------------------------------*/
/* transposed remainder of ell modulo P                                   */
/* ell has m terms, with m=deg(P)                                         */
/* res has k terms                                                        */
/*------------------------------------------------------------------------*/
void nmod_poly_trem(mp_ptr res, 
		    mp_srcptr ell, const nmod_poly_t P, const long k);

/*------------------------------------------------------------------------*/
/* transposed product: res = B.ell mod P                                  */
/* S = 1/rev(P) mod x^(m-1), with m=deg(P)                                */
/*------------------------------------------------------------------------*/
void nmod_poly_tmulmod(mp_ptr res,
		       mp_srcptr ell, const nmod_poly_t B, const nmod_poly_t P, const nmod_poly_t S);




/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* 5. to and from Newtom sums                                             */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* computes the first k Newton sums of P                                  */
/* newton must have length at least k                                     */
/*------------------------------------------------------------------------*/
void nmod_poly_to_newton_sums(mp_ptr newton, const nmod_poly_t P, long k);

/*------------------------------------------------------------------------*/
/* computes P of degree d from its first (d+1) Newton sums                */
/* d is given, since the polynomial newton may be normalized              */
/* assumes the characteristic is prime to catch error situations          */
/*------------------------------------------------------------------------*/
void nmod_poly_from_newton_sums(nmod_poly_t P, mp_srcptr newton, long d);




/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* 6. misc. calculations mod M                                            */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* returns the trace of C modulo M                                        */
/*------------------------------------------------------------------------*/
mp_limb_t nmod_poly_trace(const nmod_poly_t C, const nmod_poly_t M);

/*------------------------------------------------------------------------*/
/* given t[i]=trace(B A^i), i = 0..deg(M)-1 (where M=minpoly(A))          */
/* computes C such that B=C(A), if such a C exists                        */
/* iM = 1/M' mod M                                                        */
/*------------------------------------------------------------------------*/
void nmod_poly_convert_from_trace(nmod_poly_t C, mp_srcptr t, const nmod_poly_t M, const nmod_poly_t iM);

/*------------------------------------------------------------------------*/
/* given t[i]=trace(B A^i K^j), with                                      */
/*  i = 0..deg(M)-1 (M=minpoly(A)), j = 0..deg(N)-1 (N=minpoly(K))        */
/* computes C such that B=C(A, K), if such a C exists                     */
/* iM = 1/M' mod M and iN = 1/N' mod N                                    */
/* C is an array of size m*n (m=deg(M), n=deg(N))                         */
/* of the form C0..C_{m-1}, with Ci = coeff(C,X^i) \in Fp[y]/N(y)         */
/*------------------------------------------------------------------------*/
void nmod_poly_convert_from_trace_bi(mp_ptr C, mp_srcptr t, 
				     const nmod_poly_t M, const nmod_poly_t iM, const nmod_poly_t N, const nmod_poly_t iN);


/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* 7. geometric interpolation                                             */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* finds an element of order greater than ord modulo mod                  */
/* if no such element exists, return 0                                    */
/* if start is given, the search starts there                             */
/*------------------------------------------------------------------------*/
mp_limb_t nmod_poly_find_root_seed(long ord, long start, nmod_t mod);
mp_limb_t nmod_poly_find_root(long ord, nmod_t mod);

/*------------------------------------------------------------------------*/
/* computes the polynomial                                                */
/*     f = (x-1)(x-a)(x-a^2)...(x-a^(n-1))                                */
/*------------------------------------------------------------------------*/
void nmod_poly_eval_geom_fanin(nmod_poly_t f, mp_limb_t a, long n);

/*------------------------------------------------------------------------*/
/* Let P=(x-1)(x-a)(x-a^2)...(x-a^(n-1))                                  */
/* Evaluates P' on 1,a,...,a^(n-1)                                        */
/*------------------------------------------------------------------------*/
void nmod_poly_eval_geom_derivative_fanin(mp_ptr values, mp_limb_t a, long N, nmod_t mod);

/*------------------------------------------------------------------------*/
/*  Prepares common objects: depends only on the degree n                 */
/*  and the evaluation points. Returns                                    */
/*  1/q^{i^2} (i < 2n-1) and S=\sum_{i < 2n-1} q^{i^2} x^i                */
/*------------------------------------------------------------------------*/
void nmod_poly_eval_geom_prepare(mp_ptr inverse_powers_square_q, nmod_poly_t S, mp_limb_t q, long n);

/*------------------------------------------------------------------------*/
/* values[i] = P(q^(2i)), i = 0..n-1, where n is determined from S        */
/* inverse_powers_square_q and S are from the prepare function            */
/*------------------------------------------------------------------------*/
void nmod_poly_eval_geom_precomp(mp_ptr values, const nmod_poly_t P, mp_srcptr inverse_powers_square_q, const nmod_poly_t S);

/*------------------------------------------------------------------------*/
/* G = reverse(fanin)                                                     */
/* inverse_derivative[i] = 1/(fanin'(q^(2i)) q^(i^2), i < n               */
/*------------------------------------------------------------------------*/
void nmod_poly_interp_geom_prepare(nmod_poly_t G, mp_ptr inverse_derivative, 
				   mp_srcptr inverse_powers_square_q, mp_limb_t q, long n);

/*------------------------------------------------------------------------*/
/* finds P of degree < n such that                                        */
/* values[i] = P(q^(2i)), i = 0..n-1, where n is determined from S        */
/* all other arguments are from the prepare functions                     */
/*------------------------------------------------------------------------*/
void nmod_poly_interp_geom_precomp(nmod_poly_t P, mp_srcptr values, 
				   mp_srcptr inverse_powers_square_q, const nmod_poly_t S,
				   mp_srcptr inverse_derivative, const nmod_poly_t G);

#endif
