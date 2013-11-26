#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>

/*------------------------------------------------------------------------*/
/* computes P of degree d from its first (d+1) Newton sums                */
/* d is given, since the polynomial newton may be normalized              */
/* assumes the characteristic is prime to catch error situations          */
/*------------------------------------------------------------------------*/
void nmod_poly_from_newton_sums(nmod_poly_t P, mp_srcptr newton, long d){
  mp_limb_t n = P->mod.n;
  nmod_poly_t S;
  nmod_poly_init2(S, n, d+1);
  long i;
  mp_ptr coefS = S->coeffs;
  for (i = 0; i < d; i++)
    coefS[i] = n_negmod(newton[i+1], n);
  S->length = d;
  _nmod_poly_normalise(S);
  nmod_poly_integral(S, S);
  nmod_poly_exp_series(P, S, d+1);
  nmod_poly_reverse(P, P, d+1);
  nmod_poly_clear(S);
}
