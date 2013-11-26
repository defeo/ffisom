#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* returns the trace of C modulo M                                        */
/*------------------------------------------------------------------------*/
mp_limb_t nmod_poly_trace(const nmod_poly_t C, const nmod_poly_t M){
  nmod_poly_t dM;
  nmod_poly_init(dM, M->mod.n);
  nmod_poly_derivative(dM, M);
  nmod_poly_mulmod(dM, dM, C, M);
  mp_limb_t tmp = nmod_poly_get_coeff_ui(dM, nmod_poly_degree(M)-1);
  nmod_poly_clear(dM);
  return tmp;
}
