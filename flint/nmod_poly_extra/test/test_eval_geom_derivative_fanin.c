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
  nmod_t mod;
  nmod_init(&mod, n);

  for (i = 1; i < 200; i+=1){

    mp_ptr val = _nmod_vec_init(i);
    mp_limb_t q = nmod_poly_find_root(i, mod);

    if (opt == 1){
      nmod_poly_eval_geom_derivative_fanin(val, q, i, mod);

      sage_output_init(mod);
      printf("q = k(%lu)\n", q);
      printf("i = %lu\n", i);
      sage_output_assign_vec(val, i, "val");
      printf("P = mul([x-q^j for j in range(i)])\n");
      printf("val == [P.derivative()(q^j) for j in range(i)]\n");
    }
    else{
      double t, u;
      long j;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_eval_geom_derivative_fanin(val, q, i, mod);
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

    _nmod_vec_clear(val);
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
