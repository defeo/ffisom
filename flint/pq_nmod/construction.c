#include <stdlib.h>
#include "pq_nmod.h"
#include "nmod_poly_extra.h"
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>


void pq_nmod_init(pq_nmod_t A, const nmod_poly_t M) {
  A->degree = nmod_poly_degree(M);
  mp_limb_t n = M->mod.n;
  nmod_poly_init(A->M, n);
  nmod_poly_set(A->M, M);
  nmod_poly_init(A->iM, n);
  nmod_poly_init(A->S, n);
  nmod_poly_init(A->newton, n);
}

void _pq_nmod_init_newton(pq_nmod_t A, const nmod_t mod, 
			  mp_srcptr newton, slong k, 
			  slong degree) {
  A->degree = degree;
  nmod_poly_init(A->M, mod.n);
  nmod_poly_init(A->iM, mod.n);
  nmod_poly_init(A->S, mod.n);
  nmod_poly_init2(A->newton, mod.n, k);  
  _nmod_vec_set(A->newton->coeffs, newton, k);
  A->newton->length = k;
}

#define __NEWTON_BOUND(d, mod) (d * ((mod.n > d) + 1))

int _pq_nmod_insure_M(const pq_nmod_t A) {
  if (nmod_poly_is_zero(A->M)) {
    if (nmod_poly_is_zero(A->newton))
      return 0;
    slong 
      d = A->degree,
      k = __NEWTON_BOUND(d, A->M->mod);
    if (!A->newton->length > k) {
      printf("Exception (_pq_nmod_insure_M). Not enough Newton coefficients.\n");
      abort();
    }
    pq_nmod_struct* tmp = (pq_nmod_struct*)A;
    nmod_poly_init2(tmp->M, A->newton->mod.n, d+1);
    if (d == k)
      nmod_poly_from_newton_sums(tmp->M, A->newton->coeffs, d);
    else
      nmod_poly_minimal_polynomial_sequence(tmp->M, A->newton->coeffs, d);
  }
  return 1;
}

int _pq_nmod_insure_newton_w_length(const pq_nmod_t A, slong k) {
  if (A->newton->length < k) {
    if (nmod_poly_is_zero(A->M))
      return 0;
    pq_nmod_struct* tmp = (pq_nmod_struct*)A;
    nmod_poly_init2(tmp->newton, A->M->mod.n, k);
    // TODO: improve using trem
    nmod_poly_to_newton_sums(tmp->newton->coeffs, A->M, k);
    tmp->newton->length = k;
  }
  return 1;
}

int _pq_nmod_insure_newton(const pq_nmod_t A) {
  return _pq_nmod_insure_newton_w_length(A, A->degree);
}

int _pq_nmod_insure_iM(const pq_nmod_t A) {
  if (nmod_poly_is_zero(A->iM)) {
    if (!_pq_nmod_insure_M(A))
      return 0;
    pq_nmod_struct* tmp = (pq_nmod_struct*) A;
    nmod_poly_t dM;
    nmod_poly_init(dM, A->M->mod.n);
    nmod_poly_derivative(dM, A->M);
    if (!nmod_poly_invmod(tmp->iM, dM, A->M)) {
      printf("Exception (_pq_nmod_insure_iM). Ramified extension.\n");
      abort();
    }
    nmod_poly_clear(dM);
  }
  return 1;
}

int _pq_nmod_insure_S(const pq_nmod_t A) {
  if (nmod_poly_is_zero(A->S)) {
    if (!_pq_nmod_insure_M(A))
      return 0;
    pq_nmod_struct* tmp = (pq_nmod_struct*) A;
    nmod_poly_init(tmp->S, A->M->mod.n);
    nmod_poly_inverse_reverse_main(tmp->S, A->M);
  }
  return 1;
}

void pq_nmod_compositum(pq_nmod_t C, const pq_nmod_t A, const pq_nmod_t B) {
  slong
    m = A->degree,
    n = B->degree,
    M = m*n,
    k = __NEWTON_BOUND(M, A->M->mod);
  if (!_pq_nmod_insure_newton(A) ||
      !_pq_nmod_insure_newton(B)) {
    printf("Exception (pq_nmod_compositum). Invalid format for A or B.\n");
    abort();
  }
  nmod_poly_fit_length(C->newton, k);
  _pq_nmod_embed(C->newton->coeffs,
		 A->newton->coeffs, A->M,
		 B->newton->coeffs, B->M,
		 k);
  C->newton->length = k;
  nmod_poly_zero(C->M);
  nmod_poly_zero(C->iM);
  nmod_poly_zero(C->S);
}
