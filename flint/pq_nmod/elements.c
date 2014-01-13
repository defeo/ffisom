#include "pq_nmod.h"
#include "nmod_poly_extra.h"

inline void _pq_nmod_clear_mono(pq_nmod_elt_t x) {
  if (x->bit_field & HAS_MONO) {
    nmod_poly_clear(x->mono);
    x->bit_field &= ~HAS_MONO;
  }
}

inline void _pq_nmod_clear_dual(pq_nmod_elt_t x) {
  if (x->bit_field & HAS_DUAL) {
    nmod_poly_clear(x->dual);
    x->bit_field &= ~HAS_DUAL;
  }
}

int _pq_nmod_insure_mono(const pq_nmod_elt_t x, const pq_nmod_t A) {
  if (!(x->bit_field & HAS_MONO)) {
    if (!(x->bit_field & HAS_DUAL)) 
      return 0;
    pq_nmod_elt_struct* tmp = (pq_nmod_elt_struct*)x;
    _pq_nmod_insure_M(A);
    _pq_nmod_insure_iM(A);
    nmod_poly_init(tmp->mono, A->M->mod.n);
    nmod_poly_convert_from_trace(tmp->mono, tmp->dual->coeffs, A->M, A->iM);
    tmp->bit_field |= HAS_MONO;
  }
  return 1;
}

int _pq_nmod_insure_dual(const pq_nmod_elt_t x, const pq_nmod_t A) {
  if (!(x->bit_field & HAS_DUAL)) {
    if (!(x->bit_field & HAS_MONO)) 
      return 0;
    pq_nmod_elt_struct* tmp = (pq_nmod_elt_struct*)x;
    _pq_nmod_insure_M(A);
    _pq_nmod_insure_S(A);
    _pq_nmod_insure_trace(A);
    nmod_poly_init2(tmp->dual, A->M->mod.n, A->M->length-1);
    nmod_poly_tmulmod(tmp->dual->coeffs, A->trace, tmp->mono, A->M, A->S);
    tmp->dual->length = A->M->length-1;
    tmp->bit_field |= HAS_DUAL;
  }
  return 1;
}

