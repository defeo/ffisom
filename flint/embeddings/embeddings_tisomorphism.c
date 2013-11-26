#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------------------*/
/* input G is a linear form defined modulo R                              */
/* output F is an array of size m*n (m=deg(P), n=deg(Q))                  */
/* of the form F0..F_{m-1}, with Fi = G(S^i T^j), j=0..n-1y]/Q(y)         */
/* with S=Phi(X), T=Phi(Y), Phi = isomorphism FPxFQ -> FR                 */
/*------------------------------------------------------------------------*/

void embeddings_tisomorphism(mp_ptr F, mp_srcptr G,
			     const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR){

  long m = nmod_poly_degree(FP->P);
  long n = nmod_poly_degree(FQ->P);

  long i;
  nmod_poly_t S, X;

  nmod_t mod = FP->P->mod;
  nmod_poly_init(S, mod.n);
  nmod_poly_init(X, mod.n);

  nmod_poly_zero(X);
  nmod_poly_set_coeff_ui(X, 1, 1);
  embeddings_embed(S, X, FP, FQ, FR);

  mp_ptr tmpF = _nmod_vec_init(n);
  mp_ptr tmpG = _nmod_vec_init(m*n);

  for (i = 0; i < m*n; i++)
    tmpG[i] = G[i];

  for (i = 0; i < m; i++){
    embeddings_tembed(tmpF, tmpG, FQ, FP, FR);

    long j;
    long offset = i*n;
    for (j = 0; j < n; j++)
      F[offset+j] = tmpF[j];
    
    nmod_poly_tmulmod(tmpG, tmpG, S, FR->P, FR->SP);
  }

  _nmod_vec_clear(tmpF);
  _nmod_vec_clear(tmpG);

  nmod_poly_clear(X);
  nmod_poly_clear(S);
}
