#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <math.h>
#include "nmod_poly_extra.h"
#include "embeddings.h"
#include "sage_output.h"

/*------------------------------------------------------------------------*/
/* F is an array of size m*n (m=deg(P), n=deg(Q))                         */
/* of the form F0..F_{m-1}, with Fi = coeff(F,X^i) \in Fp[y]/Q(y)         */
/* so coeff(F, X^i Y^j) = F[i*n+j]                                        */
/*------------------------------------------------------------------------*/

void embeddings_isomorphism_2(nmod_poly_t G, mp_srcptr F, 
			      const embeddings_t FP, const embeddings_t FQ, const embeddings_t FR){

  nmod_t mod = G->mod;
  long m = nmod_poly_degree(FP->P);
  long n = nmod_poly_degree(FQ->P);

  long np = n + m - 1;
  long p = ceil(sqrt(np));
  long q = ceil((1.0*np)/p);
  long i, j, k;

  nmod_poly_t iT, iTm, T, X, tmp;
  nmod_poly_struct *TT;

  nmod_poly_init(iTm, mod.n);
  nmod_poly_init(iT, mod.n);
  nmod_poly_init(T, mod.n);
  nmod_poly_init(X, mod.n);
  nmod_poly_init(tmp, mod.n);

  nmod_poly_mat_t MT, MH, MV;
  nmod_poly_mat_init(MT, q, n, mod.n);
  nmod_poly_mat_init(MH, p, q, mod.n);
  nmod_poly_mat_init(MV, p, n, mod.n);

  TT = flint_malloc(sizeof(nmod_poly_struct) * (q+1));
  for (i = 0; i < q+1; i++)
    nmod_poly_init(TT+i, mod.n);

  nmod_poly_set_coeff_ui(X, 1, 1);
  embeddings_embed(T, X, FQ, FP, FR);
  nmod_poly_invmod(iT, T, FR->P);
  nmod_poly_powmod_ui_binexp(iTm, iT, m-1, FR->P);

  nmod_poly_zero(TT);
  nmod_poly_set_coeff_ui(TT, 0, 1);

  for (i = 1; i < q+1; i++)
    nmod_poly_mulmod(TT+i, TT+(i-1), T, FR->P);

  for (i = 0; i < q; i++)
    for (j = 0; j < n; j++){
      long jm = j*m;
      for (k = 0; k < m; k++)
  	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(MT, i, j), k, nmod_poly_get_coeff_ui(TT+i, k+jm));
    }

  for (i = 0; i < p; i++)
    for (j = 0; j < q; j++){
      long idx = i*q+j;
      long lo = FLINT_MAX(0,m-1-idx);
      long hi = FLINT_MIN(m,n+m-1-idx);
      for (k = lo; k < hi; k++)
  	nmod_poly_set_coeff_ui(nmod_poly_mat_entry(MH, i, j), k, F[k*n+k+idx-m+1]);
    }

  nmod_poly_mat_mul(MV, MH, MT);

  nmod_poly_zero(G);

  for (i = p-1; i >= 0; i--){
    nmod_poly_zero(tmp);
    for (j = 0; j < n; j++){
      long len = nmod_poly_mat_entry(MV, i, j)->length;
      mp_ptr coefs = nmod_poly_mat_entry(MV, i, j)->coeffs;
      long jm = j*m;
      for (k = 0; k < len; k++)
  	nmod_poly_set_coeff_ui(tmp, k+jm, n_addmod(nmod_poly_get_coeff_ui(tmp, k+jm), coefs[k], mod.n));
    }
    nmod_poly_rem(tmp, tmp, FR->P);
    nmod_poly_mulmod(G, G, TT+q, FR->P);
    nmod_poly_add(G, G, tmp);
  }

  nmod_poly_mulmod(G, G, iTm, FR->P);
			     

  nmod_poly_clear(tmp);
  nmod_poly_mat_clear(MT);
  nmod_poly_mat_clear(MH);
  nmod_poly_mat_clear(MV);

  for (i = 0; i < q+1; i++)
    nmod_poly_clear(TT+i);
  flint_free(TT);

  nmod_poly_clear(iTm);
  nmod_poly_clear(iT);
  nmod_poly_clear(T);
  nmod_poly_clear(X);



}
