#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"

/* ------------------------------------------------------------ */
/*  Prepares common objects: depends only on the degree n       */
/*  and the evaluation points. Returns                          */
/*  1/q^{i^2} (i < 2n-1) and S=\sum_{i < 2n-1} q^{i^2} x^i      */
/* ------------------------------------------------------------ */
void nmod_poly_eval_geom_prepare(mp_ptr inverse_powers_square_q, nmod_poly_t S, mp_limb_t q, long n){

  // computes suitable powers of q^2
  mp_limb_t power_q = 1;
  nmod_t mod = S->mod;

  nmod_poly_fit_length(S, 2*n-1);
  mp_ptr powers_square_q = S->coeffs;
  powers_square_q[0] = 1;
  inverse_powers_square_q[0] = 1;

  long i;
  for (i = 1; i < 2*n-1; i++){
    powers_square_q[i] = nmod_mul(powers_square_q[i-1], power_q, mod);
    power_q = nmod_mul(q, power_q, mod);
    powers_square_q[i] = nmod_mul(powers_square_q[i], power_q, mod);
    inverse_powers_square_q[i] = powers_square_q[i];
  }

  S->length = 2*n-1;
  _nmod_poly_normalise(S);

  _nmod_vec_invert(inverse_powers_square_q, 2*n-1, mod);
}
