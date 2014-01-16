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
  _pq_nmod_init_S(A);
  _pq_nmod_init_iM(A);
  _pq_nmod_init_newton(A);
}

void pq_nmod_init_compositum(pq_nmod_t C, const pq_nmod_t A, const pq_nmod_t B) {
  nmod_t mod = A->M->mod;
  slong
    m = A->degree,
    n = B->degree,
    d = m*n,
    k =  d * ((mod.n > d) + 1);

  C->degree = d;

  nmod_poly_init2(C->newton, mod.n, k);
  _pq_nmod_embed(C->newton->coeffs,
		 A->newton->coeffs, A->M,
		 B->newton->coeffs, B->M,
		 k);
  C->newton->length = k;

  nmod_poly_init2(C->M, mod.n, d+1);
  if (d == k)
    nmod_poly_from_newton_sums(C->M, A->newton->coeffs, d);
  else
    nmod_poly_minimal_polynomial_sequence(C->M, A->newton->coeffs, d);

  _pq_nmod_init_S(C);
  _pq_nmod_init_iM(C);
}

void _pq_nmod_init_newton(pq_nmod_t A) {
  nmod_poly_init(A->newton, A->M->mod.n);
  nmod_poly_fit_length(A->newton, A->degree);
  nmod_poly_to_newton_sums(A->newton->coeffs, A->M, A->degree);
  A->newton->length = A->degree;
}

void _pq_nmod_init_iM(pq_nmod_t A) {
  nmod_poly_t dM;
  nmod_poly_init(dM, A->M->mod.n);
  nmod_poly_derivative(dM, A->M);
  nmod_poly_init(A->iM, A->M->mod.n);
  if (!nmod_poly_invmod(A->iM, dM, A->M)) {
    printf("Exception (_pq_nmod_insure_iM). Ramified extension.\n");
    abort();
  }
  nmod_poly_clear(dM);
}

void _pq_nmod_init_S(pq_nmod_t A) {
  nmod_poly_init(A->S, A->M->mod.n);
  nmod_poly_inverse_reverse_main(A->S, A->M);
}

void pq_nmod_clear(pq_nmod_t A) {
  A->degree = 0;
  nmod_poly_clear(A->M);
  nmod_poly_clear(A->iM);
  nmod_poly_clear(A->S);
  nmod_poly_clear(A->newton);
}
