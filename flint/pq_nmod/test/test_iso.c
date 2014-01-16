#include "pq_nmod.h"
#include <flint/flint.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  slong
    i,
    p = argc > 1 ? atoi(argv[1]) : 101,
    m = argc > 2 ? atoi(argv[2]) : 2,
    n = argc > 3 ? atoi(argv[3]) : m+1;

  nmod_poly_t P, Q;
  flint_rand_t rand;
  flint_randinit(rand);
  nmod_poly_init(P, p);
  nmod_poly_init(Q, p);
  nmod_poly_randtest_monic_irreducible(P, rand, m+WORD(1));
  nmod_poly_randtest_monic_irreducible(Q, rand, n+WORD(1));
  printf("P: ");
  nmod_poly_print(P);
  printf("\nQ: ");
  nmod_poly_print(Q);

  pq_nmod_t A, B, AB;
  pq_nmod_init(A, P);
  pq_nmod_init(B, Q);
  pq_nmod_init_compositum(AB, A, B);
  printf("\nR: ");
  nmod_poly_print(AB->M);

  pq_nmod_elt_t x[n], xx[n], z;
  for (i = 0; i < n; i++) {
    pq_nmod_elt_init(x[i], p);
    pq_nmod_elt_init(xx[i], p);
    nmod_poly_randtest_monic(x[i]->mono, rand, m);
    flint_printf("\nx_%d: ", i);
    nmod_poly_print(x[i]->mono);
  }
  pq_nmod_elt_init(z, p);

  pq_nmod_iso_from_mono(z, x, A, B);
  printf("\nzT=Σ: ");
  nmod_poly_print(z->dual);

  pq_nmod_iso_to_mono(xx, A, z, AB, B);

  for (i = 0; i < n; i++) {
    if (!pq_nmod_elt_equal(x[i], xx[i], A)) {
      flint_printf("\nx_%d: ", i);
      nmod_poly_print(xx[i]->mono);
      printf("\ntest FAILED\n");
      return 0;
    }
  }
  printf("\ntest OK\n");
  return 1;
}
