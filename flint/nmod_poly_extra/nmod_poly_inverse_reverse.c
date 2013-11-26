#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* res = 1 / reverse(P) mod x^k                                           */
/*------------------------------------------------------------------------*/
void nmod_poly_inverse_reverse(nmod_poly_t res,
			       const nmod_poly_t P, const long k){
  nmod_poly_t tmp;
  if (k <= 0){
    nmod_poly_zero(res);
    return;
  }
  long m = nmod_poly_degree(P);
  nmod_poly_init2(tmp, P->mod.n, m+1);
  nmod_poly_reverse(tmp, P, m+1);
  nmod_poly_inv_series(res, tmp, k);
  nmod_poly_clear(tmp);
}


/*------------------------------------------------------------------------*/
/* res = 1 / reverse(P) mod x^{deg(P)-1}                                  */
/* (used for modular multiplication and its transpose)                    */
/*------------------------------------------------------------------------*/
void nmod_poly_inverse_reverse_main(nmod_poly_t res,
				    const nmod_poly_t P){
  nmod_poly_inverse_reverse(res, P, nmod_poly_degree(P)-1);
}


