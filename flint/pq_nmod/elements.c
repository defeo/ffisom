#include "pq_nmod.h"
#include "nmod_poly_extra.h"
#include <stdlib.h>

void pq_nmod_elt_init(pq_nmod_elt_t x, mp_limb_t n) {
  nmod_poly_init(x->mono, n);
  nmod_poly_init(x->dual, n);
}

void pq_nmod_elt_set(pq_nmod_elt_t x, const pq_nmod_elt_t y) {
  nmod_poly_set(x->mono, y->mono);
  nmod_poly_set(x->dual, y->dual);
}

void pq_nmod_elt_set_mono(pq_nmod_elt_t x, const nmod_poly_t val) {
  nmod_poly_set(x->mono, val);
  nmod_poly_zero(x->dual);
}

void pq_nmod_elt_set_dual(pq_nmod_elt_t x, const nmod_poly_t val) {
  nmod_poly_zero(x->mono);
  nmod_poly_set(x->dual, val);
}

void _pq_nmod_insure_mono(const pq_nmod_elt_t x, const pq_nmod_t A) {
  if (nmod_poly_is_zero(x->mono) &&
      !nmod_poly_is_zero(x->dual)) {
    pq_nmod_elt_struct* tmp = (pq_nmod_elt_struct*)x;
    nmod_poly_init(tmp->mono, A->M->mod.n);
    nmod_poly_convert_from_trace(tmp->mono, tmp->dual->coeffs, A->M, A->iM);
  }
}

void _pq_nmod_insure_dual(const pq_nmod_elt_t x, const pq_nmod_t A) {
  if (nmod_poly_is_zero(x->dual) &&
      !nmod_poly_is_zero(x->mono)) {
    pq_nmod_elt_struct* tmp = (pq_nmod_elt_struct*)x;
    nmod_poly_init2(tmp->dual, A->M->mod.n, A->M->length-1);
    nmod_poly_tmulmod(tmp->dual->coeffs, A->newton->coeffs, tmp->mono, A->M, A->S);
    tmp->dual->length = A->M->length-1;
  }
}

#define __DO_OP(op, res, x, y) ((x->length && y->length) ? op(res, x, y) : nmod_poly_zero(res))

void pq_nmod_add(pq_nmod_elt_t res, const pq_nmod_elt_t x, const pq_nmod_elt_t y) {
  __DO_OP(nmod_poly_add, res->mono, x->mono, y->mono);
  __DO_OP(nmod_poly_add, res->dual, x->dual, y->dual);
}

void pq_nmod_sub(pq_nmod_elt_t res, const pq_nmod_elt_t x, const pq_nmod_elt_t y) {
  __DO_OP(nmod_poly_sub, res->mono, x->mono, y->mono);
  __DO_OP(nmod_poly_sub, res->dual, x->dual, y->dual);
}

void pq_nmod_mul(pq_nmod_elt_t res, const pq_nmod_elt_t x,
		 const pq_nmod_elt_t y, const pq_nmod_t A) {
  switch (nmod_poly_is_zero(y->mono) | 
	  (nmod_poly_is_zero(y->dual) << 1) |
	  (nmod_poly_is_zero(x->mono) << 4) |
	  (nmod_poly_is_zero(x->dual) << 5)) {
    const pq_nmod_elt_struct* tmp;
  case 0x22:
    // Both have only dual -> add mono to x
    _pq_nmod_insure_mono(x, A);
  case 0x12:
  case 0x32:
    // x has mono, y has dual -> swap them
    tmp = x;
    x = y;
    y = tmp;
  case 0x21:
  case 0x23:
    // y has mono, x has dual
    nmod_poly_fit_length(res->dual, A->degree);
    nmod_poly_tmulmod(res->dual->coeffs, x->dual->coeffs, y->mono, A->M, A->S);
    nmod_poly_zero(res->mono);
    break;
  case 0x11:
  case 0x13:
  case 0x31:
  case 0x33:
    // both have mono
    nmod_poly_mulmod(res->mono, x->mono, y->mono, A->M);
    nmod_poly_zero(res->dual);
    break;
  default:
    // in any other case, result is 0
    nmod_poly_zero(res->mono);
    nmod_poly_zero(res->dual);
    break;
  }
}

int pq_nmod_inv(pq_nmod_elt_t res, const pq_nmod_elt_t x, const pq_nmod_t A) {
  _pq_nmod_insure_mono(x, A);
  if (nmod_poly_invmod(res->mono, x->mono, A->M)) {
    nmod_poly_zero(res->dual);
    return 1;
  } else {
    return 0;
  }
}

int pq_nmod_div(pq_nmod_elt_t res, const pq_nmod_elt_t x,
		const pq_nmod_elt_t y, const pq_nmod_t A) {
  const pq_nmod_elt_struct* tmp_ptr;
  pq_nmod_elt_t tmp;

  // Support for aliasing
  if (res == x) {
    pq_nmod_elt_init(tmp, A->M->mod.n);
    pq_nmod_elt_set(tmp, x);
    tmp_ptr = tmp;
  } else {
    tmp_ptr = x;
  }

  if (pq_nmod_inv(res, y, A)) {
    pq_nmod_mul(res, tmp_ptr, res, A);
    return 1;
  } else {
    return 0;
  }
}
