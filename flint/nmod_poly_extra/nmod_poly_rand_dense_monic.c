#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* random dense monic polynomial of degree len                            */
/*------------------------------------------------------------------------*/
void nmod_poly_rand_dense_monic(nmod_poly_t poly, flint_rand_t state, long len){
  nmod_poly_rand_dense(poly, state, len+1);
  nmod_poly_fit_length(poly, len+1);
  poly->coeffs[len] = 1;
  poly->length = len+1;
}
