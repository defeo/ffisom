#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* random dense polynomial with at most len terms (so degree < len)       */
/*------------------------------------------------------------------------*/
void nmod_poly_rand_dense(nmod_poly_t poly, flint_rand_t state, long len){
  int i;
  nmod_poly_fit_length(poly, len);
  for (i = 0; i < len; i++){
    poly->coeffs[i] = n_randtest(state) % poly->mod.n;
  }

  poly->length = len;
  _nmod_poly_normalise(poly);
}
