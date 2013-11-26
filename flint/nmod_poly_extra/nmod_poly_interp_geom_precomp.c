#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"


/*------------------------------------------------------------------------*/
/* finds P of degree < n such that                                        */
/* values[i] = P(q^(2i)), i = 0..n-1, where n is determined from S        */
/* all other arguments are from the prepare functions                     */
/*------------------------------------------------------------------------*/
void nmod_poly_interp_geom_precomp(nmod_poly_t P, mp_srcptr values, 
				   mp_srcptr inverse_powers_square_q, const nmod_poly_t S,
				   mp_srcptr inverse_derivative, const nmod_poly_t G){

  
  nmod_t mod = G->mod;
  long n = S->length;

  nmod_poly_t PV;
  nmod_poly_init2(PV, mod.n, n);
  long i;
  mp_ptr PVc = PV->coeffs;
  for (i = 0; i < n; i++)
    PVc[i] = values[i];
  PV->length = n;
  _nmod_poly_normalise(PV);
  nmod_poly_mullow(PV, PV, G, n);

  nmod_poly_t R;
  nmod_poly_init2(R, mod.n, n+1);
  R->coeffs[0] = 0;
  mp_ptr Rc = R->coeffs;
  for (i = 1; i <= n; i++)
    Rc[i] = nmod_mul(nmod_poly_get_coeff_ui(S, i-1), inverse_powers_square_q[i], mod);
  R->length = n+1;
  _nmod_poly_normalise(R);
  nmod_poly_mullow(R, R, S, 2*n);

  nmod_poly_fit_length(P, n);
  mp_ptr Pc = P->coeffs;
  for (i = 0; i < n; i++)
    Pc[i] = nmod_mul(nmod_poly_get_coeff_ui(R, i+n), inverse_derivative[i], mod);

  nmod_poly_clear(PV);
  nmod_poly_clear(R);
}


