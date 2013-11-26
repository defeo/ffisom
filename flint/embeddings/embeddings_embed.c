#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------------------*/
/* computes G = phi(F), where phi is the embedding FP -> FR               */
/*------------------------------------------------------------------------*/

void embeddings_embed(nmod_poly_t G, const nmod_poly_t F, 
		      const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR){

  long m = nmod_poly_degree(FP->P);
  long n = nmod_poly_degree(FQ->P);

  mp_ptr ell = _nmod_vec_init(m);
  mp_ptr traces_Q_long = _nmod_vec_init(m*n);
  mp_ptr ellstar = _nmod_vec_init(m*n);

  nmod_poly_trem(traces_Q_long, FQ->TP->coeffs, FQ->P, m*n);
  nmod_poly_tmulmod(ell, FP->TP->coeffs, F, FP->P, FP->SP);
  nmod_poly_trem(ellstar, ell, FP->P, m*n);

  long i;
  nmod_t mod = FP->P->mod;
  for (i = 0; i < m*n; i++)
    ellstar[i] = nmod_mul(ellstar[i], traces_Q_long[i], mod);

  nmod_poly_convert_from_trace(G, ellstar, FR->P, FR->iP);

  _nmod_vec_clear(traces_Q_long);
  _nmod_vec_clear(ell);
  _nmod_vec_clear(ellstar);
}
