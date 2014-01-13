#include "pq_nmod.h"
#include "nmod_poly_extra.h"
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>


void pq_nmod_init(pq_nmod_t A, slong degree) {
  A->degree = degree;
  A->bit_field = 0;
}

int _pq_nmod_insure_M(const pq_nmod_t A) {
  if (!(A->bit_field & HAS_M)) {
    if (!(A->bit_field & HAS_trace))
      return 0;
    slong 
      d = A->degree,
      n = A->degree * ((A->trace->mod.n > d) + 1);
    if (!A->trace->length > n)
      return 0;
    pq_nmod_struct* tmp = (pq_nmod_struct*)A;
    nmod_poly_init2(tmp->M, A->trace->mod.n, d+1);
    if (d == n)
      nmod_poly_from_newton_sums(tmp->M, A->trace->coeffs, d);
    else
      nmod_poly_minimal_polynomial_sequence(tmp->M, A->trace->coeffs, d);
    tmp->bit_field |= HAS_M;
  }
  return 1;
}

int _pq_nmod_insure_trace(const pq_nmod_t A) {
  if (!(A->bit_field & HAS_trace)) {
    if (!(A->bit_field & HAS_M))
      return 0;
    slong d = A->degree;
    pq_nmod_struct* tmp = (pq_nmod_struct*)A;
    nmod_poly_init2(tmp->trace, A->M->mod.n, d+1);
    nmod_poly_to_newton_sums(tmp->trace->coeffs, A->M, d);
    tmp->trace->length = d+1; 
    tmp->bit_field |= HAS_trace;
  }
  return 1;
}

int _pq_nmod_insure_iM(const pq_nmod_t A) {
  if (!(A->bit_field & HAS_iM)) {
    if (!_pq_nmod_insure_M(A))
      return 0;
    pq_nmod_struct* tmp = (pq_nmod_struct*) A;
    nmod_poly_t dM;
    nmod_poly_init(dM, A->M->mod.n);
    nmod_poly_derivative(dM, A->M);
    nmod_poly_invmod(tmp->iM, dM, A->M);
    nmod_poly_clear(dM);
    tmp->bit_field |= HAS_iM;
  }
  return 1;
}

int _pq_nmod_insure_S(const pq_nmod_t A) {
  if (!(A->bit_field & HAS_S)) {
    if (!_pq_nmod_insure_M(A))
      return 0;
    pq_nmod_struct* tmp = (pq_nmod_struct*) A;
    nmod_poly_init(tmp->S, A->M->mod.n);
    nmod_poly_inverse_reverse_main(tmp->S, A->M);
    tmp->bit_field |= HAS_S;
  }
  return 1;
}

