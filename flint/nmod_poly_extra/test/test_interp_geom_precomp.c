#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/nmod_poly.h>

#include "util.h"
#include "sage_output.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  long i;
  flint_rand_t state;
  flint_randinit(state);
  mp_limb_t n = 65537;
    
  for (i = 1; i < 100; i+=1){
    mp_ptr inverse_powers_square_q = _nmod_vec_init(2*i-1);
    mp_ptr inverse_derivative = _nmod_vec_init(i);
    nmod_poly_t S;
    nmod_poly_init(S, n);
    nmod_poly_t G;
    nmod_poly_init(G, n);
 
    mp_limb_t q = nmod_poly_find_root(2*i, S->mod);
    nmod_poly_eval_geom_prepare(inverse_powers_square_q, S, q, i);
    nmod_poly_interp_geom_prepare(G, inverse_derivative, inverse_powers_square_q, q, i);
     
    nmod_poly_t F;
    nmod_poly_init(F, n);
    nmod_poly_rand_dense(F, state, i);

    mp_ptr val = _nmod_vec_init(i);
    nmod_poly_eval_geom_precomp(val, F, inverse_powers_square_q, S);

    nmod_poly_t F2;
    nmod_poly_init(F2, n);
    
    if (opt == 1){
      nmod_poly_interp_geom_precomp(F2, val, inverse_powers_square_q, S, inverse_derivative, G);
 
      sage_output_init(S->mod);
      printf("q = k(%lu)\n", q);
      printf("i = %lu\n", i);
      sage_output_assign_poly(F, "F");
      sage_output_assign_poly(F, "F2");
      printf("F == F2\n");
    }
    else{
      double t, u, v, w;
      long j;

      // our algo
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_interp_geom_precomp(F2, val, inverse_powers_square_q, S, inverse_derivative, G);

      t = util_gettime() - t;
      
      // poly mult
      nmod_poly_t tmp1, tmp2, tmp3;
      nmod_poly_init(tmp1, n);
      nmod_poly_init(tmp2, n);
      nmod_poly_init(tmp3, n);

      nmod_poly_rand_dense(tmp1, state, i);
      nmod_poly_rand_dense(tmp2, state, i);

      u = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_mul(tmp3, tmp1, tmp2);
      u = util_gettime() - u;

      nmod_poly_clear(tmp1);
      nmod_poly_clear(tmp2);
      nmod_poly_clear(tmp3);

      // subproduct tree techniques
      mp_ptr * tree;
      mp_ptr val2 = _nmod_vec_init(i);
      mp_ptr pts = _nmod_vec_init(i);
      tree = _nmod_poly_tree_alloc(i);
      nmod_t mod = F->mod;

      _nmod_poly_tree_build(tree, pts, i, mod);

      v = util_gettime();
      for (j = 0; j < 10000; j++)
	_nmod_poly_evaluate_nmod_vec_fast_precomp(val2, F->coeffs, F->length, tree, i, mod);
      v = util_gettime() - v;

      _nmod_poly_tree_free(tree, i);
      _nmod_vec_clear(pts);
      _nmod_vec_clear(val2);

      mp_ptr val3 = _nmod_vec_init(i);
      mp_ptr pts3 = _nmod_vec_init(i);
      _nmod_vec_randtest(pts3, state, i, F->mod);

      w = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_evaluate_nmod_vec(val3, F, pts3, i);
      w = util_gettime() - w;

      _nmod_vec_clear(pts3);
      _nmod_vec_clear(val3);
      
      printf("%lu %f %f %f %f\n", i, t, u, v, w);
    }

    nmod_poly_clear(G);
    nmod_poly_clear(S);
    nmod_poly_clear(F);
    nmod_poly_clear(F2);
    _nmod_vec_clear(val);
    _nmod_vec_clear(inverse_powers_square_q);
    _nmod_vec_clear(inverse_derivative);
  }

  flint_randclear(state);
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/* if not argument is given, runs timings                     */
/* if the argument 1 is given, runs check                     */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
  int opt = 0;
  if (argc > 1)
    opt = atoi(argv[1]);
  check(opt);
  return 0;
}
