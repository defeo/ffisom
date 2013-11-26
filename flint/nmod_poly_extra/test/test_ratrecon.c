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

  for (i = 1; i < 200; i+=3){

    nmod_poly_t num1, den1, iden1, s;
    nmod_poly_init2(num1, n, i+1);
    nmod_poly_init2(den1, n, i+1);
    nmod_poly_init(s, n);
    nmod_poly_init(iden1, n);
    
    nmod_poly_randtest(num1, state, i+1);
    nmod_poly_randtest(den1, state, i+1);
    nmod_poly_set_coeff_ui(den1, 0, 2);
    
    nmod_poly_inv_series(iden1, den1, 2*i+1);
    nmod_poly_mullow(s, iden1, num1, 2*i+1);
    
    nmod_poly_t num, den;
    nmod_poly_init(num, n);
    nmod_poly_init(den, n);
    

    if (opt == 1){
      nmod_poly_ratrecon(num, den, s, i);
    
      nmod_poly_t tmp1, tmp2, s2;
      nmod_poly_init(tmp1, n);
      nmod_poly_init(tmp2, n);
      nmod_poly_init(s2, n);
      nmod_poly_mul(tmp1, num1, den);
      nmod_poly_mul(tmp2, num, den1);
      nmod_poly_sub(s2, tmp1, tmp2);
    
      sage_output_init(s2->mod);
      sage_output_assign_poly(s2, "s2");
      printf("s2 == 0\n");

      nmod_poly_clear(s2);
      nmod_poly_clear(tmp1);
      nmod_poly_clear(tmp2);
    }
    else{
      long j;
      double t, u;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_ratrecon(num, den, s, i);
      t = util_gettime() - t;
      
      u = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_mul(s, num1, den1);
      u = util_gettime() - u;
      
      printf("%lu %f %f\n", i, t, u);
    }

    nmod_poly_clear(num1);
    nmod_poly_clear(den1);
    nmod_poly_clear(iden1);
    nmod_poly_clear(num);
    nmod_poly_clear(den);
    nmod_poly_clear(s);
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
