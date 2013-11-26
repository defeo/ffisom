#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------------------*/
/* assigns all polynomials in F from the defining polynomial P            */
/*------------------------------------------------------------------------*/

void embeddings_set_parameters(embeddings_t F, const nmod_poly_t P0){
  long d = nmod_poly_degree(P0);

  nmod_poly_set(F->P, P0);

  nmod_poly_inverse_reverse_main(F->SP, F->P);

  nmod_poly_t dP;
  nmod_poly_init(dP, F->P->mod.n);
  nmod_poly_derivative(dP, F->P);
  nmod_poly_invmod(F->iP, dP, F->P);
  nmod_poly_clear(dP);
    
  mp_ptr trace_P;
  trace_P = _nmod_vec_init(d);
  nmod_poly_to_newton_sums(trace_P, F->P, d);
  nmod_poly_fit_length(F->TP, d);
  long i;
  for (i = 0; i < d; i++)
    F->TP->coeffs[i] = trace_P[i];
  _nmod_vec_clear(trace_P);

}
