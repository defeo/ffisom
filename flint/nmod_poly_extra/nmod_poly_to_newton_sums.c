#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* computes the first k Newton sums of P                                  */
/* newton must have length at least k                                     */
/*------------------------------------------------------------------------*/
void nmod_poly_to_newton_sums(mp_ptr newton, const nmod_poly_t P, long k){
  long m = nmod_poly_degree(P);
  mp_limb_t n = P->mod.n;
  nmod_poly_t tmp0, tmp1, tmp2;
  nmod_poly_init2(tmp0, n, m);
  nmod_poly_init2(tmp1, n, m);
  nmod_poly_init2(tmp2, n, k);

  nmod_poly_derivative(tmp1, P);
  nmod_poly_reverse(tmp0, tmp1, m);
  nmod_poly_inverse_reverse(tmp2, P, k);

  nmod_poly_mullow(tmp1, tmp0, tmp2, k);

  long i;
  long len1 = tmp1->length;
  mp_ptr coef1  = tmp1->coeffs;
  for (i = 0; i < len1; i++)
    newton[i] = coef1[i];
  for (i = len1; i < k; i++)
    newton[i] = 0;

  nmod_poly_clear(tmp0);
  nmod_poly_clear(tmp1);
  nmod_poly_clear(tmp2);
}
