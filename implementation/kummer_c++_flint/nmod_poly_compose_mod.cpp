#include <iostream>
#include <flint/nmod_poly.h>
#include "nmod_poly_compose_mod.h"

/*------------------------------------------------------*/
/* helper function for modular composition              */
/* computes res[i] = polys[i](arg) mod f, 0 <= i < len2 */
/*------------------------------------------------------*/
void Nmod_poly_compose_mod::nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(nmod_poly_struct * res,
										const nmod_poly_struct * polys,
										slong len2) const{
    nmod_mat_t B, C;
    mp_ptr t, h;
    slong len1;

    h = _nmod_vec_init(2*n);
    t = _nmod_vec_init(2*n);

    for (long i = 0; i < len2; i++){
      nmod_poly_init2_preinv(res + i, mod.n, mod.ninv, n);
      _nmod_poly_set_length(res + i, n);
    }

    nmod_mat_init(B, k * len2, m, mod.n);
    nmod_mat_init(C, k * len2, n, mod.n);

    /* Set rows of B to the segments of polys */
    for (long j = 0; j < len2; j++) {
        len1 = (polys + j)->length;
        for (long i = 0; i < len1 / m; i++)
            _nmod_vec_set(B->rows[i + j * k], (polys + j)->coeffs + i * m, m);
    	long i = len1 / m;
        _nmod_vec_set(B->rows[i + j * k], (polys + j)->coeffs + i * m, len1 % m);
    }

    nmod_mat_mul(C, B, A);
      
    /* Evaluate block composition using the Horner scheme */
    for (long j = 0; j < len2; j++) {
      _nmod_vec_set(t, C->rows[j * k], n);
      for (long i = n; i < 2*n; i++)
    	t[i] = 0;
      for (long i = 1; i < k; i++){
    	_nmod_poly_mul(h, C->rows[j * k + i], n, giant + (i-1)*n, n, mod);
    	_nmod_poly_add(t, h, 2*n-1, t, 2*n-1, mod);
      }
      _nmod_poly_divrem_newton_n_preinv(h, (res + j)->coeffs, t, 2*n-1, f, n+1, f_inv, len_f_inv, mod);
    }
      


    for (long i = 0; i < len2; i++)
      _nmod_poly_normalise(res + i);
    
    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(B);
    nmod_mat_clear(C);
}


/*------------------------------------------------------*/
/* setup for modular composition by arg mod poly        */
/* precomputes baby steps and giant steps               */      
/*------------------------------------------------------*/
void Nmod_poly_compose_mod::nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(mp_srcptr arg, slong len_arg,
										mp_srcptr poly, slong len_poly,
										mp_srcptr polyinv, slong len_poly_inv,
										nmod_t mod_in,
										long sz){
  n = len_poly - 1;
  f = _nmod_vec_init(len_poly);
  flint_mpn_copyi(f, poly, len_poly);

  len_f_inv = len_poly_inv;
  f_inv = _nmod_vec_init(len_f_inv);
  flint_mpn_copyi(f_inv, polyinv, len_f_inv);

  m = sz;
  k = n / m + 1;
  mod = mod_in;
  nmod_mat_init(A, m, n, mod.n);

  /* baby steps: set rows of A to powers of arg */
  A->rows[0][0] = UWORD(1);
  _nmod_vec_set(A->rows[1], arg, len_arg);
  flint_mpn_zero(A->rows[1] + len_arg, n - len_arg);
  for (long i = 2; i < m; i++)
    _nmod_poly_mulmod_preinv(A->rows[i], A->rows[i - 1], n, A->rows[1], n, f, n+1, f_inv, len_f_inv, mod);

  /* giant steps */
  giant = _nmod_vec_init(k*n);
  _nmod_poly_mulmod_preinv(giant, A->rows[m - 1], n, A->rows[1], n, f, n+1, f_inv, len_f_inv, mod);
  for (long i = 1; i < k-1; i++)
    _nmod_poly_mulmod_preinv(giant + i*n, giant + (i-1)*n, n, giant, n, f, n+1, f_inv, len_f_inv, mod);
 }

/*------------------------------------------------------*/
/* setup for modular composition by arg mod poly        */
/* precomputes baby steps and giant steps               */      
/*------------------------------------------------------*/
void Nmod_poly_compose_mod::nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(const nmod_poly_t arg,
										const nmod_poly_t poly,
										const nmod_poly_t polyinv,
										long sz){

  nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(arg->coeffs, arg->length,
						      poly->coeffs, poly->length,
						      polyinv->coeffs, polyinv->length,
						      arg->mod, sz);
}



Nmod_poly_compose_mod::~Nmod_poly_compose_mod(){
  _nmod_vec_clear(f);
  _nmod_vec_clear(f_inv);
  _nmod_vec_clear(giant);
  nmod_mat_clear(A);
}
