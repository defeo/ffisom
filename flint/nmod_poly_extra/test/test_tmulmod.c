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
  long m;
  flint_rand_t state;
  flint_randinit(state);
  mp_limb_t n = 65537;

  for (m = 1; m < 300; m += 20){
    nmod_poly_t B, P, S;

    nmod_poly_init(B, n);
    nmod_poly_init(P, n);
    nmod_poly_init(S, n);

    mp_ptr ell = _nmod_vec_init(m);
    mp_ptr ellB = _nmod_vec_init(m);

    _nmod_vec_randtest(ell, state, m, P->mod);     // len(ell) = m
    nmod_poly_rand_dense(B, state, m);     // deg(B) < m
    nmod_poly_rand_dense_monic(P, state, m); // deg(P) = m, monic

    nmod_poly_inverse_reverse_main(S, P);

    if (opt == 1){
      nmod_poly_tmulmod(ellB, ell, B, P, S);
      sage_output_init(P->mod);
      sage_output_assign_poly_from_vec(ell, m, "ell");
      sage_output_assign_poly(B, "B");
      sage_output_assign_poly_from_vec(ellB, m, "ellB");
      sage_output_assign_poly(P, "P");
      sage_output_assign_poly(S, "G");
      printf("[add([ell[i]*(B*x^j %% P)[i] for i in range(P.degree())]) for j in range(2*P.degree())] == [add([ellB[i]*(x^j %% P)[i] for i in range(P.degree())]) for j in range(2*P.degree())]\n");

    }
    else{
      double t, u;
      long j;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_tmulmod(ellB, ell, B, P, S);
      t = util_gettime() - t;

      u = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_mul(B, P, P);
      u = util_gettime() - u;

      printf("%lu %f %f\n", m, t, u);
    }
	
    _nmod_vec_clear(ell);
    _nmod_vec_clear(ellB);
    nmod_poly_clear(B);
    nmod_poly_clear(P);
    nmod_poly_clear(S);
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
