#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* computes the composed product R of P and Q                             */
/* assumes that Fp is a field and that R is irreducible                   */
/*------------------------------------------------------------------------*/

void embeddings_composed_product_small_char(nmod_poly_t R, const nmod_poly_t P, const nmod_poly_t Q){

  long m = nmod_poly_degree(P);
  long n = nmod_poly_degree(Q);
  long d = m*n;

  mp_ptr nP = _nmod_vec_init(2*d+1);
  mp_ptr nQ = _nmod_vec_init(2*d+1);
  mp_ptr nR = _nmod_vec_init(2*d+1);

  nmod_poly_to_newton_sums(nP, P, 2*d+1);
  nmod_poly_to_newton_sums(nQ, Q, 2*d+1);

  long i;
  nmod_t mod = P->mod;
  for (i = 0; i < 2*d+1; i++)
    nR[i] = nmod_mul(nP[i], nQ[i], mod);


  nmod_poly_minimal_polynomial_sequence(R, nR, d);

  _nmod_vec_clear(nP);
  _nmod_vec_clear(nQ);
  _nmod_vec_clear(nR);
}
