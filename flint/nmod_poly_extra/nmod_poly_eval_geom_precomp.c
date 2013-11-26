#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* values[i] = P(q^(2i)), i = 0..n-1, where n is determined from S        */
/* inverse_powers_square_q and S are from the prepare function            */
/*------------------------------------------------------------------------*/
void nmod_poly_eval_geom_precomp(mp_ptr values, const nmod_poly_t P, mp_srcptr inverse_powers_square_q, const nmod_poly_t S){
  
  // S has length 2n-1
  long n = (S->length+1) >> 1;

  if (n == 1){
    values[0] = nmod_poly_get_coeff_ui(P, 0);
    return;
  }

  nmod_t mod = S->mod;

  nmod_poly_t T;
  nmod_poly_init2(T, mod.n, n);
  mp_ptr cT = T->coeffs;
  long i;
  for (i = 0; i < n; i++)
    cT[i] = nmod_mul(nmod_poly_get_coeff_ui(P, n-1-i), inverse_powers_square_q[n-1-i], mod);
  T->length = n;
  _nmod_poly_normalise(T);

  // does the product
  nmod_poly_mullow(T, T, S, 2*n-1);

  // extract the output
  for (i = 0; i < n; i++)
    values[i] = nmod_mul(nmod_poly_get_coeff_ui(T, n-1+i), inverse_powers_square_q[i], mod);

  nmod_poly_clear(T);
}

