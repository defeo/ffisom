#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* computes the composed product R of P and Q                             */
/* assumes that char(Fp) > deg(P) deg(Q) and that Fp is a field           */
/*------------------------------------------------------------------------*/

void embeddings_composed_product_large_char(nmod_poly_t R, const nmod_poly_t P, const nmod_poly_t Q){

  long m = nmod_poly_degree(P);
  long n = nmod_poly_degree(Q);
  long d = m*n;

  mp_ptr nP = _nmod_vec_init(d+1);
  mp_ptr nQ = _nmod_vec_init(d+1);
  mp_ptr nR = _nmod_vec_init(d+1);

  nmod_poly_to_newton_sums(nP, P, d+1);
  nmod_poly_to_newton_sums(nQ, Q, d+1);

  long i;
  for (i = 0; i < d+1; i++)
    nR[i] = nmod_mul(nP[i], nQ[i], P->mod);

  nmod_poly_from_newton_sums(R, nR, d);

  _nmod_vec_clear(nP);
  _nmod_vec_clear(nQ);
  _nmod_vec_clear(nR);
}
