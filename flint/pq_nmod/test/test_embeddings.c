#include "pq_nmod.h"
#include <flint/flint.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  slong
    p = argc > 1 ? atoi(argv[1]) : 101,
    m = argc > 2 ? atoi(argv[2]) : 3,
    n = argc > 3 ? atoi(argv[3]) : m+1;

  nmod_poly_t P, Q;
  flint_rand_t rand;
  flint_randinit(rand);
  nmod_poly_init(P, p);
  nmod_poly_init(Q, p);
  nmod_poly_randtest_monic_irreducible(P, rand, m);
  nmod_poly_randtest_monic_irreducible(Q, rand, n);
  printf("P: ");
  nmod_poly_print(P);
  printf("\nQ: ");
  nmod_poly_print(Q);

  pq_nmod_t A, B, AB;
  pq_nmod_init(A, P);
  pq_nmod_init(B, Q);
  pq_nmod_init_compositum(AB, A, B);

  pq_nmod_elt_t x, y, z;
  pq_nmod_elt_init(x, p);
  pq_nmod_elt_init(y, p);
  pq_nmod_elt_init(z, p);
  nmod_poly_randtest_monic(P, rand, m-1);
  nmod_poly_randtest_monic(Q, rand, n-1);
  pq_nmod_elt_set_mono(x, P);
  pq_nmod_elt_set_mono(y, Q);
  printf("\nx: ");
  nmod_poly_print(x->mono);
  printf("\ny: ");
  nmod_poly_print(y->mono);

  pq_nmod_embed(z, x, A, y, B);
  printf("\nz=xy: ");
  nmod_poly_print(z->dual);

  nmod_poly_zero(P);
  nmod_poly_set_coeff_ui(P, WORD(m-1), UWORD(1));
  _pq_nmod_elt_set_dual(y, P);
  pq_nmod_project(y, A, z, AB, y, B);
  printf("\nz|(0,...,0,1): ");
  nmod_poly_print(y->mono);

  printf("\n");
  return pq_nmod_elt_equal(x, y, A);
}
