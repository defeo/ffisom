#ifndef PQ_NMOD_H
#define PQ_NMOD_H

#include <flint/nmod_poly.h>

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/


// Polynomial quotient algebras
typedef struct{
  nmod_poly_t M;
  nmod_poly_t iM;
  nmod_poly_t S;
  nmod_poly_t newton;
  slong degree;
} pq_nmod_struct;
typedef pq_nmod_struct pq_nmod_t[1];

// Elements of polynomial quotient algebras
typedef struct {
  nmod_poly_t mono;
  nmod_poly_t dual;
} pq_nmod_elt_struct;
typedef pq_nmod_elt_struct pq_nmod_elt_t[1];

#define __COEFF_PROD(res, x, y, mod, M) for ((M)--; (M) >= WORD(0); (M)--) (res)[M] = nmod_mul((x)[M], (y)[M], (mod))


/************** ALGEBRAS *****************/

void pq_nmod_init(pq_nmod_t A, const nmod_poly_t M);
void pq_nmod_init_compositum(pq_nmod_t C, const pq_nmod_t A, const pq_nmod_t B);
void _pq_nmod_init_newton(pq_nmod_t A);
void _pq_nmod_init_iM(pq_nmod_t A);
void _pq_nmod_init_S(pq_nmod_t A);
void pq_nmod_clear(pq_nmod_t A);

/************** ELEMENTS *****************/

void pq_nmod_elt_init(pq_nmod_elt_t x, mp_limb_t n);
void pq_nmod_elt_set(pq_nmod_elt_t x, const pq_nmod_elt_t y);
void pq_nmod_elt_set_mono(pq_nmod_elt_t x, const nmod_poly_t val);
void _pq_nmod_elt_set_dual(pq_nmod_elt_t x, const nmod_poly_t val);
void pq_nmod_elt_zero(pq_nmod_elt_t x);
int pq_nmod_elt_equal(const pq_nmod_elt_t x, const pq_nmod_elt_t y, const pq_nmod_t A);
void _pq_nmod_insure_mono(const pq_nmod_elt_t x, const pq_nmod_t A);
void _pq_nmod_insure_dual(const pq_nmod_elt_t x, const pq_nmod_t A);
void pq_nmod_add(pq_nmod_elt_t res, const pq_nmod_elt_t x, const pq_nmod_elt_t y);
void pq_nmod_sub(pq_nmod_elt_t res, const pq_nmod_elt_t x, const pq_nmod_elt_t y);
void pq_nmod_mul(pq_nmod_elt_t res, const pq_nmod_elt_t x,
		 const pq_nmod_elt_t y, const pq_nmod_t A);
int pq_nmod_inv(pq_nmod_elt_t res, const pq_nmod_elt_t x, const pq_nmod_t A);
int pq_nmod_div(pq_nmod_elt_t res, const pq_nmod_elt_t x,
		const pq_nmod_elt_t y, const pq_nmod_t A);

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
void _pq_nmod_project(nmod_poly_t res, const nmod_poly_t P,
		      const nmod_poly_t x,
		      mp_srcptr y, const nmod_poly_t Q);
void pq_nmod_project(pq_nmod_elt_t res, const pq_nmod_t A,
		     const pq_nmod_elt_t x, const pq_nmod_t AB,
		     const pq_nmod_elt_t y, const pq_nmod_t B);
/* Relative trace of AB/A */
void pq_nmod_trace(pq_nmod_elt_t res, const pq_nmod_t A,
		   const pq_nmod_elt_t x, const pq_nmod_t AB,
		   const pq_nmod_t B);

/************** ISOMORPHISM *****************/

/*
  res must not alias any x[i]
*/
void pq_nmod_iso_from_mono(pq_nmod_elt_t res,
			   const pq_nmod_elt_t* x, const pq_nmod_t A, 
			   const pq_nmod_t B);
/*
  x must not alias any res[i]
*/
void pq_nmod_iso_to_dual(pq_nmod_elt_t* res, const pq_nmod_t A,
			 const pq_nmod_elt_t x, const pq_nmod_t AB,
			 const pq_nmod_t B);
void pq_nmod_iso_to_mono(pq_nmod_elt_t* res, const pq_nmod_t A,
			 const pq_nmod_elt_t x, const pq_nmod_t AB,
			 const pq_nmod_t B);



#endif
