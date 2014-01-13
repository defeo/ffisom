#include "pq_nmod.h"
#include "nmod_poly_extra.h"
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>


// For private use only
// (these should go in macros, really)
static inline void __coeff_prod(mp_ptr res,
				mp_srcptr x, mp_srcptr y,
				const nmod_t mod, slong M) {
  for (M--; M >= WORD(0); M--)
    res[M] = nmod_mul(x[M], y[M], mod);
}
static inline void __zero_pad(mp_ptr res, slong start, slong end) {
  for (end--; end >= start; end--)
    res[end] = UWORD(0);
}

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
		    slong M) {
  mp_ptr tmp = _nmod_vec_init(M);
  nmod_poly_trem(tmp, y, Q, M);
  nmod_poly_trem(res, x, P, M);
  __coeff_prod(res, res, tmp, P->mod, M);
  _nmod_vec_clear(tmp);
}

/*
  Given x in A=k[X]/P and y in B=k[X]/Q, computes the image of xy in
  the compositum AB. P and Q must be monic and normalized.
 */
void pq_nmod_embed(pq_nmod_elt_t res, 
		   const pq_nmod_elt_t x, const pq_nmod_t A,
		   const pq_nmod_elt_t y, const pq_nmod_t B) {
  slong
    m = A->M->length,
    n = B->M->length,
    M = m*n;
  _pq_nmod_insure_dual(x, A);
  _pq_nmod_insure_dual(y, B);
  nmod_poly_fit_length(res->dual, M);
  _pq_nmod_embed(res->dual->coeffs, x->dual->coeffs, A->M, y->dual->coeffs, B->M, M);
  res->dual->length = M;
  _pq_nmod_clear_mono(res);
}



/************** PROJECTION *****************/

/*
  Project x onto the space (k[X]/P)y*, where y* is any element dual to
  y. x is in monomial form, y is in dual form. Result is in monomial
  form.
 */
void _pq_nmod_project(nmod_poly_t res, const nmod_poly_t x,
		      mp_srcptr y, const nmod_poly_t Q,  
		      const nmod_poly_t P) {
  slong
    m = nmod_poly_degree(P),
    n = nmod_poly_degree(Q),
    M = FLINT_MIN(x->length, m*n);
  
  mp_ptr tmp = _nmod_vec_init(M);
  nmod_poly_trem(tmp, y, Q, M);
  nmod_poly_fit_length(res, M);
  __coeff_prod(res->coeffs, x->coeffs, tmp, P->mod, M);
  res->length = M;
  _nmod_poly_normalise(res);
  _nmod_vec_clear(tmp);

  nmod_poly_rem(res, res, P);
}

void pq_nmod_project(pq_nmod_elt_t res,
		     const pq_nmod_elt_t x, const pq_nmod_t AB,
		     const pq_nmod_elt_t y, const pq_nmod_t B,
		     const pq_nmod_t A) {
  _pq_nmod_insure_mono(x, AB);
  _pq_nmod_insure_dual(y, B);
  _pq_nmod_project(res->mono, x->mono, y->mono->coeffs, B->M, A->M);
  _pq_nmod_clear_dual(res);
}
