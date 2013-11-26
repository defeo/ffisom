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

  for (i = 0; i < 1000; i+=10){

    nmod_poly_t P;
    nmod_poly_init(P, n);
    nmod_poly_rand_dense_monic(P, state, i);
    
    nmod_poly_t S0;
    nmod_poly_init(S0, n);

    if (opt == 1){
      nmod_poly_inverse_reverse(S0, P, i+10);

      sage_output_init(P->mod);
      sage_output_assign_poly(P, "P");
      sage_output_assign_poly(S0, "S0");
    
      nmod_poly_inverse_reverse_main(P, P);
      sage_output_assign_poly(P, "S1");
      printf("print (S0*P.reverse()) %% x^%lu == 1\n", i+10);
      if (i > 1)
	printf("print (S1*P.reverse()) %% x^%lu == 1\n", i-1);
    }

    else{
      double t, u;
      long j;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_inverse_reverse_main(S0, P);
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
    nmod_poly_clear(S0);
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


