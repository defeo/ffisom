#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------------------*/
/* clears all polynomials in F                                            */
/*------------------------------------------------------------------------*/

void embeddings_clear(embeddings_t F){
  nmod_poly_clear(F->P);
  nmod_poly_clear(F->SP);
  nmod_poly_clear(F->iP);
  nmod_poly_clear(F->TP);
}
