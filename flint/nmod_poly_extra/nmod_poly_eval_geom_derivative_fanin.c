#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"

/*----------------------------------------------------------------------*/
/* Let P=(x-1)(x-a)(x-a^2)...(x-a^(n-1))                                */
/* Evaluates P' on 1,a,...,a^(n-1)                                      */
/*----------------------------------------------------------------------*/
void nmod_poly_eval_geom_derivative_fanin(mp_ptr values, mp_limb_t a, long n, nmod_t mod){

  mp_ptr P = _nmod_vec_init(n);
  mp_ptr Q = _nmod_vec_init(n);
  mp_ptr powA = _nmod_vec_init(n);
  mp_ptr ipowA = _nmod_vec_init(n);
  long i;

  // powers of A
  powA[0] = 1;
  ipowA[0] = 1;
  for (i = 1; i < n; i++){
    powA[i] = nmod_mul(powA[i-1], a, mod);
    ipowA[i] = powA[i];
  }
  _nmod_vec_invert_montgomery(ipowA, n, mod);

  // small recurrence formulas give the result
  P[0] = 1;
  for (i = 1; i < n; i++){
    mp_limb_t tmp1 = nmod_mul(P[i-1], powA[i-1], mod);
    mp_limb_t tmp2 = nmod_sub(powA[i], 1, mod);
    P[i] = nmod_mul(tmp1, tmp2, mod);
  }

  Q[n-1] = 1;
  for (i = n-2; i >= 0; i--){
    mp_limb_t tmp1 = nmod_sub(powA[i], powA[n-1], mod);
    mp_limb_t tmp2 = nmod_mul(Q[i+1], tmp1, mod);
    Q[i] = nmod_mul(tmp2, ipowA[n-2-i], mod);
  }

  for (i = 0; i < n; i++)
    values[i] = nmod_mul(P[i], Q[i], mod);

  _nmod_vec_clear(ipowA);
  _nmod_vec_clear(powA);
  _nmod_vec_clear(P);
  _nmod_vec_clear(Q);
}

