#include <iostream>
#include <flint/fmpz.h>
#include <flint/fq_nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_poly.h>
#include <flint/profiler.h>
#include "nmod_poly_build_irred.h"
#include "nmod_min_poly.h"
#include "nmod_cyclotomic_poly.h"
#include "cyclotomic_ext_rth_root.h"

using namespace std;


//---------------------------------------------------
// computes F = min-poly of elt modulo <x^r+a, ctx> 
// assumes that we are in a field extension (no check is done)
// t is a bound on the degree of the minpoly we want to compute
//----------------------------------------------------
void cyclotomic_ext_min_poly_special(nmod_poly_t F, slong r, fq_nmod_t a, fq_nmod_poly_t elt, slong t, fq_nmod_ctx_t ctx){
  long len = 2*t+1;
  long k = n_sqrt(len)+1;
  long m = len / k + 1;
  long s = nmod_poly_degree(ctx->modulus);

  nmod_mat_t A, B, C;
  nmod_mat_init(A, k, r*s, F->mod.n);
  nmod_mat_init(B, r*s, m, F->mod.n);
  nmod_mat_init(C, k, m, F->mod.n);

  fq_nmod_poly_t elt_pow;
  fq_nmod_poly_init(elt_pow, ctx);
  fq_nmod_t one;
  fq_nmod_init(one, ctx);
  fq_nmod_set_si(one, 1, ctx);
  fq_nmod_poly_set_fq_nmod(elt_pow, one, ctx);

  fq_nmod_t coeff, tmp, zero;
  fq_nmod_init(coeff, ctx);
  fq_nmod_init(tmp, ctx);
  fq_nmod_init(zero, ctx);
  fq_nmod_zero(zero, ctx);

  //  timeit_t time;
  //  timeit_start(time);

  // baby steps
  // elt_pow = elt^i  
  for (long i = 0; i < k; i++){
    // put them in rows of matrix A
    for (long j = 0; j < r; j++){
      fq_nmod_poly_get_coeff(coeff, elt_pow, j, ctx);
      for (long u = 0; u < s; u++)
	nmod_mat_entry(A, i, j*s + u) = nmod_poly_get_coeff_ui(coeff, u);
    }
    // multiplication by elt; reduction done by hand, because modulus is simple
    fq_nmod_poly_mul(elt_pow, elt_pow, elt, ctx);
    for (long j = r; j <= fq_nmod_poly_degree(elt_pow, ctx); j++){
      fq_nmod_poly_get_coeff(coeff, elt_pow, j, ctx);
      fq_nmod_poly_set_coeff(elt_pow, j, zero, ctx);
      fq_nmod_mul(coeff, coeff, a, ctx);
      fq_nmod_poly_get_coeff(tmp, elt_pow, j-r, ctx);
      fq_nmod_sub(coeff, tmp, coeff, ctx);
      fq_nmod_poly_set_coeff(elt_pow, j-r, coeff, ctx);
    }
  }
  //  timeit_stop(time);
  //  cout << "baby steps: " << (double) time->wall / 1000.0 << "\n";

  // giant steps
  // pow = (last elt_pow)^i = elt^(ki)
  //  timeit_start(time);

  fq_nmod_poly_t pow;
  fq_nmod_poly_init(pow, ctx);
  fq_nmod_poly_set(pow, elt_pow, ctx);
  fq_nmod_poly_set_fq_nmod(elt_pow, one, ctx);

  for (long i = 0; i < m; i++){
    // do a transposed product of each coefficient of elt_pow with tau = [0 .. 0 1] (length s)
    // put entries in matrix B, with coefficients in reverse order w.r.t x
    for (long v = 0; v < r; v++){
      fq_nmod_poly_get_coeff(coeff, elt_pow, v, ctx);
      nmod_poly_reverse(coeff, coeff, s);
      nmod_poly_mullow(coeff, coeff, ctx->inv, s);
      for (long u = 0; u < s; u++)
	B->entries[((r-1-v)*s+u)*m + i] = nmod_poly_get_coeff_ui(coeff, u);
    }
    // multiplication by elt_pow; reduction done by hand, because modulus is simple
    fq_nmod_poly_mul(elt_pow, elt_pow, pow, ctx);
    for (long j = r; j <= fq_nmod_poly_degree(elt_pow, ctx); j++){
      fq_nmod_poly_get_coeff(coeff, elt_pow, j, ctx);
      fq_nmod_poly_set_coeff(elt_pow, j, zero, ctx);
      fq_nmod_mul(coeff, coeff, a, ctx);
      fq_nmod_poly_get_coeff(tmp, elt_pow, j-r, ctx);
      fq_nmod_sub(coeff, tmp, coeff, ctx);
      fq_nmod_poly_set_coeff(elt_pow, j-r, coeff, ctx);
    }
  }
  fq_nmod_poly_clear(pow, ctx);
  fq_nmod_poly_clear(elt_pow, ctx);
  fq_nmod_clear(one, ctx);
  fq_nmod_clear(zero, ctx);
  fq_nmod_clear(coeff, ctx);
  fq_nmod_clear(tmp, ctx);

  //  timeit_stop(time);
  //  cout << "giant steps: " << (double) time->wall / 1000.0 << "\n";

  //  timeit_start(time);
  mp_ptr seq = _nmod_vec_init(k*m);
  nmod_mat_mul(C, A, B);
  long idx = 0;
  for (long i = 0; i < m; i++)
    for (long j = 0; j < k; j++){
      seq[idx++] = nmod_mat_entry(C, j, i);
    }
  nmod_mat_clear(A);
  nmod_mat_clear(B);
  nmod_mat_clear(C);
  //  timeit_stop(time);
  //  cout << "matmul: " << (double) time->wall / 1000.0 << "\n";

  
  // minimal polynomial of the sequence, and done.
  NmodMinPoly minp;
  //  timeit_start(time);
  minp.minimal_polynomial(F, seq, t);
  //  timeit_stop(time);
  //  cout << "minpoly: " << (double) time->wall / 1000.0 << "\n";

  _nmod_vec_clear(seq);
}

