#include <flint/nmod_poly.h>
#include "nmod_poly_compose_mod.h"

static void
_nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(nmod_poly_struct * res,
						     const nmod_poly_struct * polys,
						     slong lenpolys, slong l,
						     mp_srcptr poly, slong len,
						     mp_srcptr polyinv, slong leninv,
						     nmod_t mod)
{
    nmod_mat_t A, B, C;
    mp_ptr t, h;
    slong k, n, m, len2 = l, len1;

    n = len - 1;

    m = n_sqrt(n * len2) + 1;
    k = len / m + 1;

    h = _nmod_vec_init(2*n);
    t = _nmod_vec_init(2*n + k*n);
    mp_ptr giant = t + 2*n;

    nmod_mat_init(A, m, n, mod.n);
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

    /* Set rows of A to powers of last element of polys */
    A->rows[0][0] = UWORD(1);
    _nmod_vec_set(A->rows[1], (polys + lenpolys - 1)->coeffs, (polys + lenpolys - 1)->length);
    flint_mpn_zero(A->rows[1] + (polys + lenpolys - 1)->length, n - (polys + lenpolys - 1)->length);
    for (long i = 2; i < m; i++)
        _nmod_poly_mulmod_preinv(A->rows[i], A->rows[i - 1], n, A->rows[1], n, poly, len, polyinv, leninv, mod);

    nmod_mat_mul(C, B, A);

    /* giant steps */
    _nmod_poly_mulmod_preinv(giant, A->rows[m - 1], n, A->rows[1], n, poly, len, polyinv, leninv, mod);
    for (long i = 1; i < k-1; i++)
      _nmod_poly_mulmod_preinv(giant + i*n, giant + (i-1)*n, n, giant, n, poly, len, polyinv, leninv, mod);
      
    /* Evaluate block composition using the Horner scheme */
    for (long j = 0; j < len2; j++) {
        _nmod_vec_set(t, C->rows[j * k], n);
	for (long i = n; i < 2*n; i++)
	  t[i] = 0;
        for (long i = 1; i < k; i++){
	  _nmod_poly_mul(h, C->rows[j * k + i], n, giant + (i-1)*n, n, mod);
	  _nmod_poly_add(t, h, 2*n-1, t, 2*n-1, mod);
        }
	_nmod_poly_divrem_newton_n_preinv(h, (res + j)->coeffs, t, 2*n-1, poly, len, polyinv, leninv, mod);
    }

    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}


/*------------------------------------------------------------------------*/
/* multiple modular composition                                           */
/* res[i] = polys[i](polys[len1-1]) mod poly, i=0..n-1                    */
/* polyinv = 1 /rev(poly) mod x^n (or n-1?)                               */
/* requires n < len1                                                      */
/* the entries of res need to be uninitialised                            */
/*------------------------------------------------------------------------*/
void Nmod_poly_compose_mod::nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(nmod_poly_struct * res,
										const nmod_poly_struct * polys,
										slong len1, slong n,
										const nmod_poly_t poly,
										const nmod_poly_t polyinv){
  slong len2 = poly->length;
  if (len2 < 10 || n < 2){
    nmod_poly_compose_mod_brent_kung_vec_preinv(res, polys, len1, n, poly, polyinv);
    return;
  }

  for (long i = 0; i < len1; i++){
    slong len3 = (polys + i)->length;
    if (len3 >= len2){
      flint_printf
	("Exception (nmod_poly_compose_mod_brent_kung_vec_preinv_precomp)."
	 "The degree of the first polynomial must be smaller than that of the "
	 " modulus\n");
      abort();
        }
    }

  if (n > len1){
    flint_printf
      ("Exception (nmod_poly_compose_mod_brent_kung_vec_preinv_precomp)."
       "n is larger than the length of polys\n");
    abort();
  }

  for (long i = 0; i < n; i++){
    nmod_poly_init2_preinv(res + i, poly->mod.n, poly->mod.ninv, len2 - 1);
    _nmod_poly_set_length(res + i, len2 - 1);
  }
  
  _nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(res, polys, len1, n, poly->coeffs, len2, polyinv->coeffs, polyinv->length, poly->mod);

  for (long i = 0; i < n; i++)
    _nmod_poly_normalise(res + i);
}

