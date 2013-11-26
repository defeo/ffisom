#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------------------*/
/* initializes all polynomials in F                                       */
/*------------------------------------------------------------------------*/

void embeddings_init(embeddings_t F, mp_limb_t n){
  nmod_poly_init(F->P, n);
  nmod_poly_init(F->SP, n);
  nmod_poly_init(F->iP, n);
  nmod_poly_init(F->TP, n);
}
