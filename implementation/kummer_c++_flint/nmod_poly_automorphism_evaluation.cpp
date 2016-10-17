#include <iostream>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_compose_mod.h"
#include "nmod_poly_automorphism_evaluation.h"

using namespace std;


void Nmod_poly_automorphism_evaluation::_automorphism_evaluation_compose(mp_ptr res, 
									 mp_srcptr a, slong len_a,
									 mp_srcptr g, 
									 mp_srcptr f, slong len_f,
									 mp_srcptr f_inv, slong len_f_inv, nmod_t mod){


    nmod_mat_t A, B, C;
    mp_ptr t, h;
    slong i, n, m;

    n = len_f - 1;
    m = 0.5*n_sqrt(len_a) + 1;
    //    m = len_a;
    //    m = len_a/2;
    long k = (len_a-1) / m + 1;

    nmod_mat_init(A, m, n, mod.n);
    nmod_mat_init(B, k, m, mod.n);
    nmod_mat_init(C, k, n, mod.n);

    h = _nmod_vec_init(n);
    t = _nmod_vec_init(n);

    /* Set rows of B to the segments of a */
    for (i = 0; i < k-1; i++)
        _nmod_vec_set(B->rows[i], a + i*m, m);
    _nmod_vec_set(B->rows[i], a + i*m, (len_a - (k-1)*m));

    /* Set rows of A to frobeniuses of g */
    flint_mpn_copyi(A->rows[0], g, n);
    for (i = 1; i < m; i++)
      _nmod_poly_powmod_ui_binexp_preinv (A->rows[i], A->rows[i-1], mod.n, 
					  f, len_f, f_inv, len_f_inv, mod);

    nmod_mat_mul(C, B, A);

    // get x^{p^m} mod f
    _nmod_poly_powmod_x_ui_preinv(h, mod.n, f, len_f, f_inv, len_f_inv, mod);
    for (i = 1; i < m; i++){
      _nmod_poly_powmod_ui_binexp_preinv (t, h, mod.n, f, len_f, f_inv, len_f_inv, mod);
      _nmod_vec_set(h, t, n);
    }

    /* Evaluate block composition using the Horner scheme */
    
    Nmod_poly_compose_mod compose;
    compose.nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(h, n, f, len_f, f_inv, len_f_inv, mod, 2*n_sqrt(n)+1);

    nmod_poly_t input;
    nmod_poly_init2_preinv(input, mod.n, mod.ninv, n);
    nmod_poly_t output;
    _nmod_vec_set(input->coeffs, C->rows[k - 1], n);
    nmod_poly_init2_preinv(output, mod.n, mod.ninv, n);
    input->length = n;
    _nmod_poly_normalise(input);
    
    for (i = k- 2; i >= 0; i--){
      compose.nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(output, input, 1);
      _nmod_poly_add(input->coeffs, output->coeffs, n, C->rows[i], n, mod);
      input->length = n;
      _nmod_poly_normalise(input);
    }

    _nmod_vec_set(res, input->coeffs, n);
    nmod_poly_clear(output);
    nmod_poly_clear(input);

    // _nmod_vec_set(res, C->rows[k - 1], n);
    // for (i = k- 2; i >= 0; i--){
    //   _nmod_poly_compose_mod_brent_kung_preinv(t, res, n, h, f, len_f, f_inv, len_f_inv, mod);
    // }


    _nmod_vec_clear(h);
    _nmod_vec_clear(t);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
}


void Nmod_poly_automorphism_evaluation::compose(nmod_poly_t res,
						const nmod_poly_t A, const nmod_poly_t g, const nmod_poly_t f, const nmod_poly_t f_inv){

    slong len_A = A->length;
    slong len_g = g->length;
    slong len_f = f->length;
    slong len = len_f - 1;

    mp_ptr ptr2;

    if (len_f == 0) {
        flint_printf("Exception (Nmod_poly_automorphism_evaluation::compose). Division by zero.\n");
        abort();
    }

    if (len_A == 0 || len_f == 1) {
        nmod_poly_zero(res);
        return;
    }

    if (len_A == 1){
	nmod_poly_scalar_mul_nmod(res, g, A->coeffs[0]);
        return;
    }

    if (res == A || res == g || res == f || res == f_inv){
        flint_printf("Exception (Nmod_poly_automorphism_evaluation::compose). Aliasing not supported.\n");
        abort();
    }

    ptr2 = _nmod_vec_init(len);

    if (len_g <= len){
        flint_mpn_copyi(ptr2, g->coeffs, len_g);
        flint_mpn_zero(ptr2 + len_g, len - len_g);
    }
    else {
        _nmod_poly_rem(ptr2, g->coeffs, len_g, f->coeffs, len_f, res->mod);
    }

    nmod_poly_fit_length(res, len);
    _automorphism_evaluation_compose(res->coeffs,
				     A->coeffs, len_A, 
				     ptr2, 
				     f->coeffs, len_f,
				     f_inv->coeffs, f_inv->length, res->mod);
    res->length = len;
    _nmod_poly_normalise(res);
    _nmod_vec_clear(ptr2);
}


/*
 * Asserts that degree(A) >= 0.
 */
void Nmod_poly_automorphism_evaluation::compose_naive(nmod_poly_t res,
						      const nmod_poly_t A, const nmod_poly_t g, const nmod_poly_t f, const nmod_poly_t f_inv){

  nmod_poly_t g_loc, tmp;
  nmod_poly_init(g_loc, f->mod.n);
  nmod_poly_set(g_loc, g);
  nmod_poly_init(tmp, f->mod.n);

  {
    nmod_poly_scalar_mul_nmod(tmp, g_loc, nmod_poly_get_coeff_ui(A, 0));
    nmod_poly_set(res, tmp);
  }
  for (slong i = 1; i <= nmod_poly_degree(A); i++) {
    nmod_poly_powmod_ui_binexp_preinv(g_loc, g_loc, f->mod.n, f, f_inv);
    nmod_poly_scalar_mul_nmod(tmp, g_loc, nmod_poly_get_coeff_ui(A, i));
    nmod_poly_add(res, res, tmp);
  }

  nmod_poly_clear(tmp);
  nmod_poly_clear(g_loc);
}
