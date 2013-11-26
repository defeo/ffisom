#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* minimal polynomial of degree at most d of the sequence val             */
/* val must have 2d+1 (or more) entries                                   */
/* (in reality, 2d would be enough but the code doesn't support it)       */
/*------------------------------------------------------------------------*/
void nmod_poly_minimal_polynomial_sequence(nmod_poly_t res, mp_srcptr val, long d){
  nmod_t mod = res->mod;
  nmod_poly_t num, den, s;
  nmod_poly_init(num, mod.n);
  nmod_poly_init(den, mod.n);

  nmod_poly_init2(s, mod.n, 2*d+1);
  mp_ptr coefS = s->coeffs;
  long i;
  for (i = 0; i < 2*d+1; i++)
    coefS[i] = val[i];
  s->length = 2*d+1;
  _nmod_poly_normalise(s);

  nmod_poly_ratrecon(num, den, s, d);
  nmod_poly_reverse(res, den, den->length);
  nmod_poly_make_monic(res, res);

  nmod_poly_clear(s);
  nmod_poly_clear(num);
  nmod_poly_clear(den);
}
