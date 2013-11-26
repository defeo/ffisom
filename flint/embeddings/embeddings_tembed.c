#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------------------*/
/* transpose version of the following:                                    */
/* computes G = phi(F), where phi is the embedding FP -> FR               */
/* Equivalently, computes ell(S^i mod R), for i < m*n, where S=phi(X)     */
/*------------------------------------------------------------------------*/
void embeddings_tembed(mp_ptr ellstar, mp_srcptr ell, 
		       const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR){

  long m = nmod_poly_degree(FP->P);
  long n = nmod_poly_degree(FQ->P);

  mp_ptr traces_Q_long = _nmod_vec_init(m*n);
  nmod_poly_trem(traces_Q_long, FQ->TP->coeffs, FQ->P, m*n);

  nmod_t mod = FP->P->mod;
  nmod_poly_t pV;
  nmod_poly_init2(pV, mod.n, m*n);

  nmod_poly_convert_from_trace(pV, ell, FR->P, FR->iP);

  mp_ptr pVc = pV->coeffs;
  long i;
  long d = pV->length;
  for (i = 0; i < d; i++)
    pVc[i] = nmod_mul(pVc[i], traces_Q_long[i], mod);

  nmod_poly_rem(pV, pV, FP->P);

  nmod_poly_tmulmod(ellstar, FP->TP->coeffs, pV, FP->P, FP->SP);

  nmod_poly_clear(pV);
  _nmod_vec_clear(traces_Q_long);
}
