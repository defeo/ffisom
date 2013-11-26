#include <flint/nmod_poly.h>

/*------------------------------------------------------------------------*/
/* transposed product of ell by poly2                                     */
/* input: len(ell) = k+m, deg(poly2) <= m                                 */
/* output: len(res) = k                                                   */
/*------------------------------------------------------------------------*/
void nmod_poly_tmul(mp_ptr res, 
		    mp_srcptr ell, const nmod_poly_t poly2, const long m, const long k){
  nmod_poly_t poly1, tmp1, tmp2;
  mp_limb_t n = poly2->mod.n;
  long i;

  nmod_poly_init2(tmp1, n, m+1);
  nmod_poly_init2(poly1, n, k+m);
  nmod_poly_init2(tmp2, n, k+m);

  mp_ptr coef1 = poly1->coeffs;
  for (i = 0; i < k+m; i++)
    coef1[i] = ell[i];
  poly1->length = k+m;
  _nmod_poly_normalise(poly1);

  nmod_poly_reverse(tmp1, poly2, m+1);
  nmod_poly_mullow(tmp2, poly1, tmp1, k+m);

  long len2 = (tmp2->length) - m;
  if (len2 < 0)
    len2 = 0;
  mp_ptr coef2 = tmp2->coeffs;
  for (i = 0; i < len2; i++)
    res[i] = coef2[i+m];
  for (i = len2; i < k; i++)
    res[i] = 0;

  nmod_poly_clear(poly1);
  nmod_poly_clear(tmp1);
  nmod_poly_clear(tmp2);
}
