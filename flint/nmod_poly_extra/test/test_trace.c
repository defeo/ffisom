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

  for (i = 1; i < 1000; i+=10){

    mp_limb_t n = 0;
    while (n < i)
      n = n_randtest_prime(state, 0);

    nmod_poly_t P;
    mp_ptr roots;
    nmod_poly_t A;
    
    nmod_poly_init(P, n);
    roots = _nmod_vec_init(i);
    _nmod_vec_randtest(roots, state, i, P->mod);
    nmod_poly_product_roots_nmod_vec(P, roots, i);
    nmod_poly_init(A, n);
    nmod_poly_rand_dense(A, state, i);

    if (opt == 1){
      mp_limb_t tr1 = nmod_poly_trace(A, P);
      mp_limb_t tr2 = 0;
      long j;
      for (j = 0; j < i; j++)
      	tr2 = nmod_add(tr2,  nmod_poly_evaluate_nmod(A, roots[j]), P->mod);

      sage_output_init(P->mod);
      sage_output_assign_poly(P, "P");
      printf("tr1 = %lu\n", tr1);
      printf("tr2 = %lu\n", tr2);
      printf("tr1 == tr2\n");
    }
    else{
      double t, u;
      long j;

      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_trace(A, P);
      t = util_gettime() - t;

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

      printf("%lu %f %f\n", i, t, u);
    }

    nmod_poly_clear(P);
    nmod_poly_clear(A);
    _nmod_vec_clear(roots);

  }
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


