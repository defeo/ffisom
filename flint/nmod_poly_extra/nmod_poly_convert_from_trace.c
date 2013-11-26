#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* given t[i]=trace(B A^i), i = 0..deg(M)-1 (where M=minpoly(A))          */
/* computes C such that B=C(A), if such a C exists                        */
/* iM = 1/M' mod M                                                        */
/*------------------------------------------------------------------------*/
void nmod_poly_convert_from_trace(nmod_poly_t C, mp_srcptr t, const nmod_poly_t M, const nmod_poly_t iM){
  long m = nmod_poly_degree(M);
  mp_limb_t n = M->mod.n;

  nmod_poly_t revM, N, PT;
  nmod_poly_init2(revM, n, m+1);
  nmod_poly_init2(N, n, m);
  nmod_poly_init2(PT, n, m);

  nmod_poly_reverse(revM, M, m+1);

  long j;
  for (j = 0; j < m; j++)
    PT->coeffs[j] = t[j];
  PT->length = m;
  _nmod_poly_normalise(PT);

  nmod_poly_mullow(N, revM, PT, m);
  nmod_poly_reverse(N, N, m);

  nmod_poly_mulmod(C, N, iM, M);

  nmod_poly_clear(revM);
  nmod_poly_clear(PT);
  nmod_poly_clear(N);
}
