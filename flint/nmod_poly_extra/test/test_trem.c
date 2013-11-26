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
  
  for (m = 1; m < 400; m += 13){
    nmod_poly_t P;
    mp_ptr ell, res;

    ell = _nmod_vec_init(m);
    res = _nmod_vec_init(2*m+5); 
    nmod_poly_init(P, n);
    
    _nmod_vec_randtest(ell, state, m, P->mod); // m terms in ell
    nmod_poly_rand_dense_monic(P, state, m); // deg(P) = m, monic

    if (opt == 1){
      nmod_poly_trem(res, ell, P, 2*m+5);
      sage_output_init(P->mod);
      sage_output_assign_poly_from_vec(ell, m, "ell");
      sage_output_assign_poly_from_vec(res, 2*m+5, "res");
      sage_output_assign_poly(P, "P");
      printf("k = %lu\n", 2*m+5);
      printf("print all([add([res[j+i]*P[j] for j in range(P.degree()+1)])==0 for i in range(k-P.degree())])\n");
    }
    else{
      double t, u;
      long j;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_trem(res, ell, P, 2*m);
      t = util_gettime() - t;

      nmod_poly_t tmp1, tmp2, tmp3;
      nmod_poly_init(tmp1, n);
      nmod_poly_init(tmp2, n);
      nmod_poly_init(tmp3, n);

      nmod_poly_rand_dense(tmp1, state, m);
      nmod_poly_rand_dense(tmp2, state, m);

      u = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_mul(tmp3, tmp1, tmp2);
      u = util_gettime() - u;

      nmod_poly_clear(tmp1);
      nmod_poly_clear(tmp2);
      nmod_poly_clear(tmp3);

      printf("%lu %f %f\n", m, t, u);
    }
	
    _nmod_vec_clear(ell);
    _nmod_vec_clear(res);
    nmod_poly_clear(P);
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
