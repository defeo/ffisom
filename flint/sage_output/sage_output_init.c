#include <stdio.h>
#include <flint/nmod_poly.h>

/*---------------------------------------------------------------*/
/* assign k=Z/pZ and U=k[x]                                      */
/*---------------------------------------------------------------*/
void sage_output_init(nmod_t p){
  printf("p = %lu\n", p.n);
  printf("k = Integers(p)\n");
  printf("U.<x> = PolynomialRing(k)\n");
}
