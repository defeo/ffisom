#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/nmod_poly.h>

#include "util.h"
#include "sage_output.h"
#include "nmod_vec_extra.h"
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
    long degP = i;
    long degQ = i+1;

    nmod_poly_t P;
    nmod_poly_init(P, n);
    nmod_poly_rand_dense_monic(P, state, degP);

    nmod_poly_t iP, dP;
    nmod_poly_init(iP, n);
    nmod_poly_init(dP, n);
    nmod_poly_derivative(dP, P);
    nmod_poly_invmod(iP, dP, P);
    nmod_poly_clear(dP);


    nmod_poly_t Q;
    nmod_poly_init(Q, n);
    nmod_poly_rand_dense_monic(Q, state, degQ);

    nmod_poly_t iQ, dQ;
    nmod_poly_init(iQ, n);
    nmod_poly_init(dQ, n);
    nmod_poly_derivative(dQ, Q);
    nmod_poly_invmod(iQ, dQ, Q);
    nmod_poly_clear(dQ);


    mp_ptr A_trace, A_poly;
    A_trace = _nmod_vec_init(degP*degQ);
    A_poly = _nmod_vec_init(degP*degQ);

    nmod_vec_rand_dense(A_trace, state, degP*degQ, P->mod);


    if (opt == 1){
      nmod_poly_convert_from_trace_bi(A_poly, A_trace, P, iP, Q, iQ);

      sage_output_init(P->mod);
      sage_output_assign_vec(A_trace, degP*degQ, "Atrace");
      sage_output_assign_vec(A_poly, degP*degQ, "Apoly");
      sage_output_assign_poly(P, "P");
      sage_output_assign_poly(Q, "Q");

      printf("M.<X,Y,Z> = PolynomialRing(GF(k.cardinality()), 3, order='lex')\n");
      printf("I = Ideal([P(X), Q(Y), Z-X*Y])\n");
      printf("GB = I.groebner_basis()\n");
      printf("quo = M.quotient(I)\n");
      printf("R = GB[2](0,0,x)\n");
      printf("qU = U.quo(R)\n");
      printf("S = qU( (quo(X)).lift()(0,0,x))\n");
      printf("T = qU( (quo(Y)).lift()(0,0,x))\n");
      printf("ApolyXY = add([Apoly[i*Q.degree()+j]*X^i*Y^j for i in range(P.degree()) for j in range(Q.degree())])\n");
      printf("A = qU( (quo(ApolyXY)).lift()(0,0,x))\n");
      printf("Atrace == [(A*S^i*T^j).trace() for i in range(P.degree()) for j in range(Q.degree())]\n");
    }

    _nmod_vec_clear(A_trace);
    _nmod_vec_clear(A_poly);

    nmod_poly_clear(P);
    nmod_poly_clear(iP);
    nmod_poly_clear(Q);
    nmod_poly_clear(iQ);

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
