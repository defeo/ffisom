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
    nmod_poly_init(P, n);
    nmod_poly_rand_dense_monic(P, state, i);

    nmod_poly_t S;
    nmod_poly_init(S, n);
    nmod_poly_inverse_reverse_main(S, P);

    nmod_poly_t iP, dP;
    nmod_poly_init(iP, n);
    nmod_poly_init(dP, n);
    nmod_poly_derivative(dP, P);
    nmod_poly_invmod(iP, dP, P);
    nmod_poly_clear(dP);

    mp_ptr trace_vec;
    trace_vec = _nmod_vec_init(i);
    nmod_poly_to_newton_sums(trace_vec, P, i);

    nmod_poly_t A;
    nmod_poly_init(A, n);
    nmod_poly_rand_dense(A, state, i);

    mp_ptr A_trace;
    A_trace = _nmod_vec_init(i);
    nmod_poly_tmulmod(A_trace, trace_vec, A, P, S);

    nmod_poly_t A2;
    nmod_poly_init(A2, n);

    if (opt == 1){
      nmod_poly_convert_from_trace(A2, A_trace, P, iP);

      sage_output_init(P->mod);
      sage_output_assign_vec(A_trace, i, "Atrace");
      sage_output_assign_poly(P, "P");
      sage_output_assign_poly(A, "A");
      sage_output_assign_poly(A2, "A2");
      printf("A == A2\n");
    }
    else{
      double t, u;
      long j;

      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_convert_from_trace(A2, A_trace, P, iP);
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
    nmod_poly_clear(iP);
    nmod_poly_clear(S);
    nmod_poly_clear(A);
    nmod_poly_clear(A2);
    _nmod_vec_clear(trace_vec);
    _nmod_vec_clear(A_trace);
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


