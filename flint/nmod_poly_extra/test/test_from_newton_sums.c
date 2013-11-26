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

  for (i = 1; i < 17; i+=1){
    nmod_poly_t P;
    mp_ptr newtonP, roots;

    nmod_poly_init(P, n);
    roots = _nmod_vec_init(i);
    _nmod_vec_randtest(roots, state, i, P->mod);
    nmod_poly_product_roots_nmod_vec(P, roots, i);

    newtonP = _nmod_vec_init(i+1);
    long k, ell;
    for (k = 0; k < i+1; k++){
      mp_limb_t s = 0;
      for (ell = 0; ell < i; ell++)
	s = n_addmod(s, n_powmod2(roots[ell], k, n), n);
      newtonP[k] = s;
    }

    if (opt == 1){
      nmod_poly_from_newton_sums(P, newtonP, i);
      
      sage_output_init(P->mod);
      sage_output_assign_poly(P, "P");
      sage_output_assign_vec(newtonP, i, "newtonP");
      sage_output_assign_vec(roots, i, "roots");
      printf("P == mul([x-r for r in roots])\n");
    }
    else{
      double t, u;
      long j;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_from_newton_sums(P, newtonP, i);
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
    _nmod_vec_clear(newtonP);
    _nmod_vec_clear(roots);
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


