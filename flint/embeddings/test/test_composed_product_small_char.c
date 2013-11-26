#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/nmod_poly.h>

#include "util.h"
#include "sage_output.h"
#include "nmod_poly_extra.h"
#include "embeddings.h"

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  long i;
  flint_rand_t state;
  flint_randinit(state);
  mp_limb_t n = 5;

  for (i = 1; i < 400; i+=3){

    nmod_poly_t P, Q, R;

    nmod_poly_init(P, n);
    do
      nmod_poly_rand_dense_monic(P, state, i);
    while(! nmod_poly_is_irreducible(P));

    nmod_poly_init(Q, n);
    do
      nmod_poly_rand_dense_monic(Q, state, i+1);
    while(! nmod_poly_is_irreducible(Q));

    nmod_poly_init(R, n);

    if (opt == 1){
      sage_output_init(P->mod);
      embeddings_composed_product_small_char(R, P, Q);
      sage_output_assign_poly(R, "R");
      sage_output_assign_poly(P, "P");
      sage_output_assign_poly(Q, "Q");
      printf("K.<X,Y,Z> = PolynomialRing(GF(P.parent().base_ring().cardinality()), 3, order='lex')\n");
      printf("I = [P(X), Q(Y), Z-X*Y]\n");
      printf("GB = Ideal(I).groebner_basis()\n");
      printf("GB[2] == R(Z)\n");
    }
    else{
      double t, u;
      long j;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	embeddings_composed_product_large_char(R, P, Q);
      t = util_gettime() - t;

      nmod_poly_t tmp1, tmp2, tmp3;
      nmod_poly_init(tmp1, n);
      nmod_poly_init(tmp2, n);
      nmod_poly_init(tmp3, n);

      nmod_poly_rand_dense(tmp1, state, i*(i+1));
      nmod_poly_rand_dense(tmp2, state, i*(i+1));

      u = util_gettime();
      for (j = 0; j < 10000; j++)
	nmod_poly_mul(tmp3, tmp1, tmp2);
      u = util_gettime() - u;

      nmod_poly_clear(tmp1);
      nmod_poly_clear(tmp2);
      nmod_poly_clear(tmp3);

      printf("%lu %f %f\n", i, t, u);
    }

    nmod_poly_clear(R);
    nmod_poly_clear(P);
    nmod_poly_clear(Q);
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
