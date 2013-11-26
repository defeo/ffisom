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
    nmod_poly_t b;
    mp_ptr a, c;

    nmod_poly_init(b, n);
    nmod_poly_rand_dense(b, state, i+1); // deg(b) = (m=i)

    a = _nmod_vec_init(2*i); // len(a) = (k+m) = (i+i) 
    _nmod_vec_randtest(a, state, 2*i, b->mod);

    c = _nmod_vec_init(i); 
 
    if (opt == 1){
      nmod_poly_tmul(c, a, b, i, i);
      sage_output_init(b->mod);
      printf("i=%ld\n", i);
      sage_output_assign_poly_from_vec(a, 2*i, "a");
      sage_output_assign_poly(b, "b");
      sage_output_assign_poly_from_vec(c, i, "c");
      printf("a.parent()([add([b[u]*a[u+j] for u in range(i+1)]) for j in range(i)]) == c\n");
    }
    else{
      double t, u;
      long j;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_tmul(c, a, b, i, i);
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
	
    nmod_poly_clear(b);
    _nmod_vec_clear(a);
    _nmod_vec_clear(c);
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
