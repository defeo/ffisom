#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------------------*/
/* F is an array of size m*n (m=deg(P), n=deg(Q))                         */
/* of the form F0..F_{m-1}, with Fi = coeff(F,X^i) \in Fp[y]/Q(y)         */
/*------------------------------------------------------------------------*/

void embeddings_isomorphism(nmod_poly_t G, mp_srcptr F, 
			    const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR){

  long m = nmod_poly_degree(FP->P);
  long n = nmod_poly_degree(FQ->P);

  long i;
  nmod_poly_t tmpF, tmpG, S, X;

  nmod_t mod = FP->P->mod;
  nmod_poly_init(tmpF, mod.n);
  nmod_poly_init(tmpG, mod.n);
  nmod_poly_init(S, mod.n);
  nmod_poly_init(X, mod.n);

  nmod_poly_zero(G);
  nmod_poly_zero(X);
  nmod_poly_set_coeff_ui(X, 1, 1);
  embeddings_embed(S, X, FP, FQ, FR);
  
  for (i = m-1; i >= 0; i--){
    nmod_poly_fit_length(tmpF, n);
    long j;
    long offset = i*n;
    for (j = 0; j < n; j++)
      tmpF->coeffs[j] = F[offset+j];
    tmpF->length = n;
    _nmod_poly_normalise(tmpF);
    
    embeddings_embed(tmpG, tmpF, FQ, FP, FR);
    nmod_poly_mulmod(G, G, S, FR->P);
    nmod_poly_add(G, G, tmpG);
  }

  nmod_poly_clear(tmpF);
  nmod_poly_clear(tmpG);
  nmod_poly_clear(X);
  nmod_poly_clear(S);
}
