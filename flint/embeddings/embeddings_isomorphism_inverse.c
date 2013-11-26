#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------------------*/
/* F is an array of size m*n (m=deg(P), n=deg(Q))                         */
/* of the form F0..F_{m-1}, with Fi = coeff(F,X^i) \in Fp[y]/Q(y)         */
/*------------------------------------------------------------------------*/

void embeddings_isomorphism_inverse(mp_ptr F, const nmod_poly_t G, 
				    const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR){

  long m = nmod_poly_degree(FP->P);
  long n = nmod_poly_degree(FQ->P);

  mp_ptr ell = _nmod_vec_init(m*n);
  nmod_poly_tmulmod(ell, FR->TP->coeffs, G, FR->P, FR->SP);

  mp_ptr Ftrace = _nmod_vec_init(m*n);
  embeddings_tisomorphism(Ftrace, ell, FP, FQ, FR);

  nmod_poly_convert_from_trace_bi(F, Ftrace, FP->P, FP->iP, FQ->P, FQ->iP);

  _nmod_vec_clear(Ftrace);
  _nmod_vec_clear(ell);
}
