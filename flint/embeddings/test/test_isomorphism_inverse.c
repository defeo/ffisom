#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/nmod_poly.h>

#include "util.h"
#include "sage_output.h"
#include "nmod_vec_extra.h"
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

  mp_limb_t n = 65537;

  for (i = 1; i*(i+1) < 500; i+=1){

    long degP = i, degQ = i+1;

    embeddings_t FP;
    embeddings_init(FP, n);
    nmod_poly_t P;
    nmod_poly_init(P, n);
    do
      nmod_poly_rand_dense_monic(P, state, degP);
    while (! nmod_poly_is_irreducible(P));
    embeddings_set_parameters(FP, P);

    embeddings_t FQ;
    embeddings_init(FQ, n);
    nmod_poly_t Q;
    nmod_poly_init(Q, n);
    do
      nmod_poly_rand_dense_monic(Q, state, degQ);
    while (! nmod_poly_is_irreducible(Q));
    embeddings_set_parameters(FQ, Q);
    
    embeddings_t FR;
    embeddings_init(FR, n);
    nmod_poly_t R;
    nmod_poly_init(R, n);
    embeddings_composed_product_small_char(R, P, Q);
    embeddings_set_parameters(FR, R);

    nmod_poly_clear(P);
    nmod_poly_clear(Q);
    nmod_poly_clear(R);

    mp_ptr Fptr = _nmod_vec_init(degP*degQ);
    mp_ptr Fptr2 = _nmod_vec_init(degP*degQ);
    nmod_vec_rand_dense(Fptr, state, degP*degQ, FP->P->mod);

    nmod_poly_t G;
    nmod_poly_init(G, n);

    embeddings_isomorphism(G, Fptr, FP, FQ, FR);
    

    if (opt == 1){
      embeddings_isomorphism_inverse(Fptr2, G, FP, FQ, FR);

      sage_output_init(P->mod);
      sage_output_assign_poly(FP->P, "P");
      sage_output_assign_poly(FQ->P, "Q");
      sage_output_assign_poly(FR->P, "R");
      sage_output_assign_vec(Fptr, degP*degQ, "F");
      sage_output_assign_vec(Fptr2, degP*degQ, "F2");
      sage_output_assign_poly(G, "G");

      printf("M.<X,Y,Z> = PolynomialRing(GF(k.cardinality()), 3, order='lex')\n");
      printf("I = Ideal([P(X), Q(Y), Z-X*Y])\n");
      printf("GB = I.groebner_basis()\n");
      printf("quo = M.quotient(I)\n");
      printf("PF = add([U([F[i+j*Q.degree()] for i in range(Q.degree())])(Y)*X^j for j in range(P.degree())])\n");
      printf("quo(PF).lift() == G(Z) and F == F2\n");
    }
    else{
      double t, u;
      long j;

      t = util_gettime();
      for (j = 0; j < 1000; j++)
	embeddings_isomorphism_inverse(Fptr2, G, FP, FQ, FR);
      t = util_gettime() - t;

      nmod_poly_t tmp1, tmp2, tmp3;
      nmod_poly_init(tmp1, n);
      nmod_poly_init(tmp2, n);
      nmod_poly_init(tmp3, n);

      nmod_poly_rand_dense(tmp1, state, degP*degQ);
      nmod_poly_rand_dense(tmp2, state, degP*degQ);

      u = util_gettime();
      for (j = 0; j < 1000; j++)
	nmod_poly_mul(tmp3, tmp1, tmp2);
      u = util_gettime() - u;

      nmod_poly_clear(tmp1);
      nmod_poly_clear(tmp2);
      nmod_poly_clear(tmp3);

      printf("%lu %f %f\n", i, t, u);
    }

    embeddings_clear(FP);
    embeddings_clear(FQ);
    embeddings_clear(FR);

    _nmod_vec_clear(Fptr2);
    _nmod_vec_clear(Fptr);
    nmod_poly_clear(G);

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


