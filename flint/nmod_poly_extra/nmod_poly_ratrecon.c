#include <flint/nmod_poly.h>

/*------------------------------------------------------------------------*/
/* rational reconstruction                                                */
/* s must have degree <= 2d                                               */
/* output: deg(num) <= d, deg(den) <= d                                   */
/* BUGBUGBUGBUGBUG: sigma                                                 */
/*------------------------------------------------------------------------*/
void nmod_poly_ratrecon(nmod_poly_t num, nmod_poly_t den, const nmod_poly_t s, long d){
  long i, ell = 2*d+1;
  nmod_t smod = s->mod;
  long sn = smod.n, sninv = smod.ninv;

  nmod_poly_t G, m0, m1, m2, m3, tnum, tden;

  nmod_poly_init2(G, sn, ell+1);
  nmod_poly_init2(m0, sn, ell+1);
  nmod_poly_init2(m1, sn, ell+1);
  nmod_poly_init2(m2, sn, ell+1);
  nmod_poly_init2(m3, sn, ell+1);
  nmod_poly_init2(tnum, sn, ell+1);
  nmod_poly_init2(tden, sn, ell+1);

  for (i = 0; i < ell; i++)
    nmod_poly_set_coeff_ui(G, i, 0);
  nmod_poly_set_coeff_ui(G, ell, 1);

  long lA, lB, lM[4];
  mp_ptr Mat[4] = {tden->coeffs, m1->coeffs, m2->coeffs, m3->coeffs};

  _nmod_poly_hgcd(Mat, lM, m0->coeffs, &lA, tnum->coeffs, &lB, G->coeffs, G->length, s->coeffs, s->length, smod);

  tden->length = lM[0];
  _nmod_poly_normalise(tden);
  m1->length = lM[1];
  _nmod_poly_normalise(m1);
  m2->length = lM[2];
  _nmod_poly_normalise(m2);
  m3->length = lM[3];
  _nmod_poly_normalise(m3);
  m0->length = lA;
  _nmod_poly_normalise(m0);
  tnum->length = lB;
  _nmod_poly_normalise(tnum);

  long det = n_submod(n_mulmod2_preinv(nmod_poly_get_coeff_ui(tden, 0), nmod_poly_get_coeff_ui(m3, 0), sn, sninv),
		      n_mulmod2_preinv(nmod_poly_get_coeff_ui(m1, 0), nmod_poly_get_coeff_ui(m2, 0), sn, sninv),
		      sn);

  /* shouldn't we use sigma? (the output of hgcd)  */
  /* if (sigma < 0)                                */
  /*   det = n_submod(0, det, sn);                 */
  
  nmod_poly_scalar_mul_nmod(den, tden, det);
  nmod_poly_set(num, tnum);

  nmod_poly_clear(G);
  nmod_poly_clear(tnum);
  nmod_poly_clear(tden);
  nmod_poly_clear(m0);
  nmod_poly_clear(m1);
  nmod_poly_clear(m2);
  nmod_poly_clear(m3);
}
