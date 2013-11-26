#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* transposed product: res = B.ell mod P                                  */
/* S = 1/rev(P) mod x^(m-1), with m=deg(P)                                */
/*------------------------------------------------------------------------*/
void nmod_poly_tmulmod(mp_ptr res,
		       mp_srcptr ell, const nmod_poly_t B, const nmod_poly_t P, const nmod_poly_t S){

  long m = nmod_poly_degree(P);
  long i;

  // init
  nmod_poly_t tmp2;
  long n = P->mod.n;
  nmod_poly_init2(tmp2, n, m-1);

  // first tmul; could be a short mul?
  mp_ptr res_loc = _nmod_vec_init(2*m-1);
  for (i = 0; i < m; i++)
    res_loc[i] = ell[i];
  for (i = m; i < 2*m-1; i++)
    res_loc[i] = 0;
  nmod_poly_tmul(tmp2->coeffs, res_loc, P, m, m-1);
  tmp2->length = m-1;
  _nmod_poly_normalise(tmp2);

  // short mul
  nmod_poly_mullow(tmp2, S, tmp2, m-1);

  // set res_loc = ell cat -tmp2
  for (i = 0; i < m; i++)
    res_loc[i] = ell[i];
  long len2 = tmp2->length;
  mp_ptr coef2 = tmp2->coeffs;
  for (i = 0; i < len2; i++)
    res_loc[i+m] = n_negmod(coef2[i], n);
  for (i = len2; i < m-1; i++)
    res_loc[i+m] = 0;

  nmod_poly_tmul(res, res_loc, B, m-1, m);

  nmod_poly_clear(tmp2);
  _nmod_vec_clear(res_loc);
}
