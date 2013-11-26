#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------------------*/
/* given t[i]=trace(B A^i K^j), with                                      */
/*  i = 0..deg(M)-1 (M=minpoly(A)), j = 0..deg(N)-1 (N=minpoly(K))        */
/* computes C such that B=C(A, K), if such a C exists                     */
/* iM = 1/M' mod M and iN = 1/N' mod N                                    */
/* C is an array of size m*n (m=deg(M), n=deg(N))                         */
/* of the form C0..C_{m-1}, with Ci = coeff(C,X^i) \in Fp[y]/N(y)         */
/*------------------------------------------------------------------------*/
void nmod_poly_convert_from_trace_bi(mp_ptr C, mp_srcptr t, 
				     const nmod_poly_t M, const nmod_poly_t iM, const nmod_poly_t N, const nmod_poly_t iN){

  long m = nmod_poly_degree(M);
  long n = nmod_poly_degree(N);
  nmod_t mod = M->mod;

  nmod_poly_t revM, revN;
  nmod_poly_init(revM, mod.n);
  nmod_poly_reverse(revM, M, m+1);
  nmod_poly_init(revN, mod.n);
  nmod_poly_reverse(revN, N, n+1);

  nmod_poly_t tmpF, tmpC;
  nmod_poly_init(tmpC, mod.n);
  nmod_poly_init(tmpF, mod.n);

  long i, j;
  for (i = 0; i < m; i++){
    nmod_poly_fit_length(tmpF, n);
    long offset = i*n;
    for (j = 0; j < n; j++)
      tmpF->coeffs[j] = t[offset+j];
    tmpF->length = n;
    _nmod_poly_normalise(tmpF);
    
    nmod_poly_mullow(tmpC, revN, tmpF, n);
    nmod_poly_reverse(tmpC, tmpC, n);
    nmod_poly_mulmod(tmpC, tmpC, iN, N);
    
    for (j = 0; j < tmpC->length; j++)
      C[offset+j] = tmpC->coeffs[j];
    for (j = tmpC->length; j < n; j++)
      C[offset+j] = 0;
  }

  for (i = 0; i < n; i++){
    nmod_poly_fit_length(tmpF, m);
    for (j = 0; j < m; j++)
      tmpF->coeffs[j] = C[i+j*n];
    tmpF->length = m;
    _nmod_poly_normalise(tmpF);
    
    nmod_poly_mullow(tmpC, revM, tmpF, m);
    nmod_poly_reverse(tmpC, tmpC, m);
    nmod_poly_mulmod(tmpC, tmpC, iM, M);
    
    for (j = 0; j < tmpC->length; j++)
      C[i+j*n] = tmpC->coeffs[j];
    for (j = tmpC->length; j < m; j++)
      C[i+j*n] = 0;
  }

  nmod_poly_clear(tmpF);
  nmod_poly_clear(tmpC);
  nmod_poly_clear(revM);
  nmod_poly_clear(revN);

}
