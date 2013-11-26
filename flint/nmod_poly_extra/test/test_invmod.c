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

  for (i = 1; i < 200; i+=20){
    nmod_poly_t A, B, C;
    nmod_poly_init(A, n);
    nmod_poly_init(B, n);
    nmod_poly_init(C, n);

    nmod_poly_rand_dense_monic(A, state, i); 
    nmod_poly_rand_dense(B, state, i); 

    if (opt == 1){
      nmod_poly_invmod(C, B, A);
      sage_output_init(A->mod);
      sage_output_assign_poly(A, "A");
      sage_output_assign_poly(B, "B");
      sage_output_assign_poly(C, "C");
      printf("C*B %% A == 1\n");
    }
    else{
      double t, u;
      long j;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_invmod(C, B, A);
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
	
    nmod_poly_clear(A);
    nmod_poly_clear(B);
    nmod_poly_clear(C);
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
