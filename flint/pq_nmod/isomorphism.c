#include "pq_nmod.h"
#include "nmod_poly_extra.h"
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>

/*
  res must not alias any x[i]
*/
void pq_nmod_iso_from_mono(pq_nmod_elt_t res,
			   const pq_nmod_elt_t* x, const pq_nmod_t A, 
			   const pq_nmod_t B) {
  slong
    m = nmod_poly_degree(A->M),
    n = nmod_poly_degree(B->M),
    M = m*n;
  int res_is_zero = 1;
  _pq_nmod_set_newton_length(B, m+n-1);
  mp_ptr tmp1 = _nmod_vec_init(M);
  nmod_poly_fit_length(res->dual, M);
  _nmod_vec_zero(res->dual->coeffs, M);
  for (n--; n >= 0; n--) {
    _pq_nmod_insure_dual(x[n], A);
    if (!nmod_poly_is_zero(x[n]->dual)) {
      _pq_nmod_embed(tmp1, x[n]->dual->coeffs, A->M, 
		     B->newton->coeffs + n, B->M, M);
      _nmod_vec_add(res->dual->coeffs, res->dual->coeffs, tmp1, M, A->M->mod);
      res_is_zero = 0;
    }
  }
  res->dual->length = res_is_zero ? 0 : M;
  nmod_poly_zero(res->mono);
}

/*
  x must not alias any res[i]
*/
void pq_nmod_iso_to_dual(pq_nmod_elt_t* res,
			 const pq_nmod_elt_t x, const pq_nmod_t AB,
			 const pq_nmod_t B, const pq_nmod_t A) {
  slong
    m = nmod_poly_degree(A->M),
    n = nmod_poly_degree(B->M);
  _pq_nmod_set_newton_length(B, m+n-1);
  _pq_nmod_insure_mono(x, AB);
  for (n--; n >= 0; n--) {
    _pq_nmod_project(res[n]->mono, x->mono, 
		     B->newton->coeffs + n, B->M, A->M);
    nmod_poly_zero(res[n]->dual);
  }
}


void pq_nmod_iso_to_mono(pq_nmod_elt_t* res,
			 const pq_nmod_elt_t x, const pq_nmod_t AB,
			 const pq_nmod_t B, const pq_nmod_t A) {
}

void pq_nmod_iso_from_BSGS(pq_nmod_elt_t res,
			   const pq_nmod_elt_t* x, const pq_nmod_t A, 
			   const pq_nmod_t B) {
}
void pq_nmod_iso_to_BSGS(pq_nmod_elt_t* res,
			 const pq_nmod_elt_t x, const pq_nmod_t AB,
			 const pq_nmod_t B, const pq_nmod_t A) {
}


