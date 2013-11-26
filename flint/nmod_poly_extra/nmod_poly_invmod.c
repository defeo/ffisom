#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* sets iA to 1/A mod M                                                   */
/* there is no check that the inverse exists!                             */
/*------------------------------------------------------------------------*/
void nmod_poly_invmod(nmod_poly_t iA, const nmod_poly_t A, const nmod_poly_t P){
  mp_limb_t n = P->mod.n;
  nmod_poly_t dummy1, dummy2;
  nmod_poly_init(dummy1, n);
  nmod_poly_init(dummy2, n);
  nmod_poly_xgcd(dummy1, dummy2, iA, P, A);
  nmod_poly_clear(dummy1);
  nmod_poly_clear(dummy2);
}
