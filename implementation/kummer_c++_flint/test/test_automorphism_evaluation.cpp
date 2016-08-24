#include <iostream>
#include <flint/nmod_poly.h>
#include "nmod_poly_automorphism_evaluation.h"
#include <flint/profiler.h>

using namespace std;

/*------------------------------------------------------------------------*/
/* random dense polynomial with at most len terms (so degree < len)       */
/*------------------------------------------------------------------------*/
void nmod_poly_rand_dense(nmod_poly_t poly, flint_rand_t state, long len){
  int i;
  nmod_poly_fit_length(poly, len);
  for (i = 0; i < len; i++){
    poly->coeffs[i] = n_randtest(state) % poly->mod.n;
  }

  poly->length = len;
  _nmod_poly_normalise(poly);
}


/*------------------------------------------------------------------------*/
/* random dense monic polynomial of degree len                            */
/*------------------------------------------------------------------------*/
void nmod_poly_rand_dense_monic(nmod_poly_t poly, flint_rand_t state, long len){
  nmod_poly_rand_dense(poly, state, len+1);
  nmod_poly_fit_length(poly, len+1);
  poly->coeffs[len] = 1;
  poly->length = len+1;
}

/*------------------------------------------------------------------------*/
/* checks automorphism evaluation                                         */
/*------------------------------------------------------------------------*/
void test_automorphism_evaluate(slong degree) {
  cout << degree << " ";
  mp_limb_t p = 9001;
  
  flint_rand_t state;
  flint_randinit(state);

  nmod_poly_t A, f, finv, g;
  nmod_poly_init(A, p);
  nmod_poly_rand_dense(A, state, degree);
  nmod_poly_init(f, p);
  nmod_poly_rand_dense_monic(f, state, degree);
  nmod_poly_init(finv, p);
  nmod_poly_reverse(finv, f, f->length);
  nmod_poly_inv_series(finv, finv, f->length);
  nmod_poly_init(g, p);
  nmod_poly_rand_dense(g, state, degree);


  Nmod_poly_automorphism_evaluation eval;
  nmod_poly_t res1, res2;
  nmod_poly_init(res1, p);
  nmod_poly_init(res2, p);

  timeit_t time;

  timeit_start(time);
  eval.compose(res1, A, g, f, finv);
  timeit_stop(time);
  cout << (double) time->wall / 1000.0 << " ";

  timeit_start(time);
  eval.compose_naive(res2, A, g, f, finv);
  timeit_stop(time);
  cout << (double) time->wall / 1000.0;

  cout << endl;

  nmod_poly_sub(res1, res1, res2);
  if (!nmod_poly_is_zero(res1))
    cout << "oops\n";

  
  nmod_poly_clear(res1);
  nmod_poly_clear(res2);
  
  flint_randclear(state);
  nmod_poly_clear(A);
  nmod_poly_clear(f);
  nmod_poly_clear(g);
}


int main() {
  for (slong i = 307; i < 307+1; i++) {
    test_automorphism_evaluate(i);
  }
  
  return 0;
}

