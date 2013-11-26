#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------------------*/
/* given G, computes F such that G = Phi(F),                              */
/*  where Phi is the isomorphism FPxFQ -> FR                              */
/* F is an array of size m*n (m=deg(P), n=deg(Q))                         */
/* of the form F0..F_{m-1}, with Fi = coeff(F,X^i) \in Fp[y]/Q(y)         */
/*------------------------------------------------------------------------*/

void embeddings_embed_inverse(nmod_poly_t F, const nmod_poly_t G, 
			      const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR){

  long m = nmod_poly_degree(FP->P);
  long n = nmod_poly_degree(FQ->P);

  mp_ptr ell = _nmod_vec_init(m*n);
  nmod_poly_tmulmod(ell, FR->TP->coeffs, G, FR->P, FR->SP);

  mp_ptr v = _nmod_vec_init(m);
  embeddings_tembed(v, ell, FP, FQ, FR);

  nmod_poly_convert_from_trace(F, v, FP->P, FP->iP);
  mp_limb_t invn = n_invmod(n, FP->P->mod.n);
  
  nmod_poly_scalar_mul_nmod(F, F, invn);

  _nmod_vec_clear(ell);
  _nmod_vec_clear(v);
}
