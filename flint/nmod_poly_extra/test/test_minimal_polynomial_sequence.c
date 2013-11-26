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
  long i, j;
  flint_rand_t state;
  flint_randinit(state);
  mp_limb_t n = 65537;
  
  for (i = 2; i < 200; i+=3){

    // the minimal polynomial of this sequence is x^2-x-1
    mp_ptr(val) = _nmod_vec_init(2*i+1);
    val[0] = 0;
    val[1] = 1;
    for (j = 2 ; j < 2*i+1; j++)
      val[j] = n_addmod(val[j-1], val[j-2], n);
    
    nmod_poly_t minpoly;
    nmod_poly_init(minpoly, n);
   
    if (opt == 1){
      nmod_poly_minimal_polynomial_sequence(minpoly, val, i);
      sage_output_init(minpoly->mod);
      sage_output_assign_poly(minpoly, "minpoly");
      printf("minpoly == (x^2-x-1)\n");
    }
    else{
      double t, u;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_minimal_polynomial_sequence(minpoly, val, i);
      t = util_gettime() - t;

      nmod_poly_t a, b, c;
      nmod_poly_init(a, n);
      nmod_poly_init(b, n);
      nmod_poly_init(c, n);

      nmod_poly_rand_dense(a, state, i);
      nmod_poly_rand_dense(b, state, i);
      
      u = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_mul(c, a, b);
      u = util_gettime() - u;
      
      printf("%lu %f %f\n", i, t, u);

      nmod_poly_clear(a);
      nmod_poly_clear(b);
      nmod_poly_clear(c);
    }
    
    nmod_poly_clear(minpoly);
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
