#ifndef PQ_NMOD_H
#define PQ_NMOD_H

#include <flint/nmod_poly.h>

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/


// Polynomial quotient algebras
#define HAS_M 1
#define HAS_iM 2
#define HAS_S 4
#define HAS_trace 8
typedef struct{
  nmod_poly_t M;
  nmod_poly_t iM;
  nmod_poly_t S;
  nmod_poly_t newton;
  slong degree;
} pq_nmod_struct;
typedef pq_nmod_struct pq_nmod_t[1];

// Elements of polynomial quotient algebras
#define HAS_MONO 1
#define HAS_DUAL 2
typedef struct {
  nmod_poly_t mono;
  nmod_poly_t dual;
  short bit_field;
} pq_nmod_elt_struct;
typedef pq_nmod_elt_struct pq_nmod_elt_t[1];


/************** ALGEBRAS *****************/

void pq_nmod_init(pq_nmod_t A, const nmod_poly_t M);
void _pq_nmod_init_newton(pq_nmod_t A, const nmod_t mod, 
			  mp_srcptr newton, slong k, 
			  slong degree);
int _pq_nmod_insure_M(const pq_nmod_t A);
int _pq_nmod_insure_newton_w_length(const pq_nmod_t A, slong k);
int _pq_nmod_insure_newton(const pq_nmod_t A);
int _pq_nmod_insure_iM(const pq_nmod_t A);
int _pq_nmod_insure_S(const pq_nmod_t A);
void pq_nmod_compositum(pq_nmod_t C, const pq_nmod_t A, const pq_nmod_t B);

/************** ELEMENTS *****************/

inline void _pq_nmod_clear_mono(pq_nmod_elt_t x);
inline void _pq_nmod_clear_dual(pq_nmod_elt_t x);
int _pq_nmod_insure_mono(const pq_nmod_elt_t x, const pq_nmod_t A);
int _pq_nmod_insure_dual(const pq_nmod_elt_t x, const pq_nmod_t A);

/************** EMBEDDING *****************/

/*
  Given x in A=k[X]/P and y in B=k[X]/Q, computes the image of xy in
  the compositum AB. x, y and the result are in dual form.
  
  M is the product of deg P and deg Q. res must have enough memory to
  hold M limbs.

  res must not equal to P->coeffs. You would be crazy to do so.
 */
void _pq_nmod_embed(mp_ptr res,
		    mp_srcptr x, const nmod_poly_t P, 
		    mp_srcptr y, const nmod_poly_t Q,
		    slong M);
/*
  Given x in A=k[X]/P and y in B=k[X]/Q, computes the image of xy in
  the compositum AB. P and Q must be monic and normalized.
 */
void pq_nmod_embed(pq_nmod_elt_t res, 
		   const pq_nmod_elt_t x, const pq_nmod_t P,
		   const pq_nmod_elt_t y, const pq_nmod_t Q);


/************** PROJECTION *****************/

/*
  Project x onto the space (k[X]/P)y*, where y* is any element dual to
  y. x is in monomial form, y is in dual form. Result is in monomial
  form.
 */
void _pq_nmod_project(nmod_poly_t res, const nmod_poly_t x,
		      mp_srcptr y, const nmod_poly_t Q,  
		      const nmod_poly_t P);
void pq_nmod_project(pq_nmod_elt_t res,
		     const pq_nmod_elt_t x, const pq_nmod_t R,
		     const pq_nmod_elt_t y, const pq_nmod_t Q,
		     const pq_nmod_t P);

#endif
