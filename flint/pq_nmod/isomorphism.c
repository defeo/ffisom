#include "pq_nmod.h"
#include "nmod_poly_extra.h"
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>

/*
  res must not alias any x[i]
*/
void pq_nmod_iso_from_mono(pq_nmod_elt_t res,
			    pq_nmod_elt_t* x, const pq_nmod_t A, 
			   const pq_nmod_t B) {
  slong
    k,
    m = nmod_poly_degree(A->M),
    n = nmod_poly_degree(B->M),
    M = m*n;
  int res_is_zero = 1;

  mp_ptr y = _nmod_vec_init(M+n-1);
  nmod_poly_trem(y, B->newton->coeffs, B->M, M+n-1);
  mp_ptr tmp = _nmod_vec_init(M);
  nmod_poly_fit_length(res->dual, M);
  _nmod_vec_zero(res->dual->coeffs, M);
  for (n--; n >= 0; n--) {
    _pq_nmod_insure_dual(x[n], A);
    if (!nmod_poly_is_zero(x[n]->dual)) {
      nmod_poly_trem(tmp, x[n]->dual->coeffs, A->M, M);
      k = M;
      __COEFF_PROD(tmp, tmp, y+n, A->M->mod, k);
      _nmod_vec_add(res->dual->coeffs, res->dual->coeffs, tmp, M, A->M->mod);
      res_is_zero = 0;
    }
  }
  _nmod_vec_clear(y);
  _nmod_vec_clear(tmp);
  res->dual->length = res_is_zero ? 0 : M;
  nmod_poly_zero(res->mono);
}

/*
  x must not alias any res[i]
*/
void pq_nmod_iso_to_dual(pq_nmod_elt_t* res, const pq_nmod_t A,
			 const pq_nmod_elt_t x, const pq_nmod_t AB,
			 const pq_nmod_t B) {
  slong
    k,
    n = nmod_poly_degree(B->M),
    M = nmod_poly_degree(AB->M);
  mp_ptr y = _nmod_vec_init(M+n-1);
  nmod_poly_trem(y, B->newton->coeffs, B->M, M+n-1);
  _pq_nmod_insure_mono(x, AB);
  for (n--; n >= 0; n--) {
    nmod_poly_fit_length(res[n]->mono, M);
    k = FLINT_MIN(x->mono->length, M);
    res[n]->mono->length = k;
    __COEFF_PROD(res[n]->mono->coeffs, x->mono->coeffs, y+n, A->M->mod, k);
    _nmod_poly_normalise(res[n]->mono);
    nmod_poly_rem(res[n]->mono, res[n]->mono, A->M);
    nmod_poly_zero(res[n]->dual);
  }
  _nmod_vec_clear(y);
}


/*
  x must not alias any res[i]
*/
void pq_nmod_iso_to_mono(pq_nmod_elt_t* res, const pq_nmod_t A,
			 const pq_nmod_elt_t x, const pq_nmod_t AB,
			 const pq_nmod_t B) {
  slong
    k,
    n = nmod_poly_degree(B->M),
    M = nmod_poly_degree(AB->M);
  mp_ptr y = _nmod_vec_init(M+n-1);
  _nmod_vec_zero(y, n);
  _pq_nmod_insure_mono(x, AB);
  for (n--; n >= 0; n--) {
    y[n] = WORD(1);
    nmod_poly_trem(y, y, B->M, M);
    nmod_poly_fit_length(res[n]->mono, M);
    k = FLINT_MIN(x->mono->length, M);
    res[n]->mono->length = k;
    __COEFF_PROD(res[n]->mono->coeffs, x->mono->coeffs, y, A->M->mod, k);
    _nmod_poly_normalise(res[n]->mono);
    nmod_poly_rem(res[n]->mono, res[n]->mono, A->M);
    nmod_poly_zero(res[n]->dual);
    y[n] = WORD(0);
  }
  _nmod_vec_clear(y);
}
