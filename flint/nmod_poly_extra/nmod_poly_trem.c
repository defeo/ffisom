#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* transposed remainder of ell modulo P                                   */
/* ell has m terms, with m=deg(P)                                         */
/* res has k terms                                                        */
/*------------------------------------------------------------------------*/
void nmod_poly_trem(mp_ptr res, 
		    mp_srcptr ell, const nmod_poly_t P, const long k){

  long m = nmod_poly_degree(P);
  long n = P->mod.n;
  long i;

  if (k <= m){
    for (i = 0; i < k; i++)
      res[i] = ell[i];
    return;
  }

  nmod_poly_t tmp1, tmp3;
  nmod_poly_init2(tmp1, n, k-m);
  nmod_poly_init2(tmp3, n, k);
  nmod_poly_inverse_reverse(tmp1, P, k-m);

  // this is rather stupid: we pad ell with zeros and do a tmul
  // this should become a short mul, somehow.

  // we take care of aliasing
  if (res == P->coeffs){
    nmod_poly_t P2;
    nmod_poly_init(P2, n);
    nmod_poly_set(P2, P);
    for (i = 0; i < m; i++)
      res[i] = ell[i];
    for (i = m; i < k; i++)
      res[i] = 0;
    nmod_poly_tmul(tmp3->coeffs, res, P2, m, k-m);
    nmod_poly_clear(P2);
  }
  else{
    for (i = 0; i < m; i++)
      res[i] = ell[i];
    for (i = m; i < k; i++)
      res[i] = 0;
    nmod_poly_tmul(tmp3->coeffs, res, P, m, k-m);
  }
  tmp3->length = k-m;
  _nmod_poly_normalise(tmp3);

  nmod_poly_mullow(tmp3, tmp1, tmp3, k-m);

  mp_ptr coef3 = tmp3->coeffs;
  long len3 = tmp3->length;
  for (i = 0; i < len3; i++)
    res[i+m] = n_negmod(coef3[i], n);

  for (i = len3; i < k-m; i++)
    res[i+m] = 0;

  nmod_poly_clear(tmp1);
  nmod_poly_clear(tmp3);
}
