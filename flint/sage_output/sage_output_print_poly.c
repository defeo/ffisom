#include <stdio.h>
#include <flint/nmod_poly.h>

/*---------------------------------------------------------------*/
/* prints P                                                      */
/*---------------------------------------------------------------*/
void sage_output_print_poly(nmod_poly_t P){
  long i;
  long m = nmod_poly_degree(P);

  if (m == -1){
    printf("0");
    return;
  }

  for (i = 0; i < m; i++)
    printf("%lu*x^%ld+", P->coeffs[i], i);
  printf("%lu*x^%ld", P->coeffs[m], m);
}
