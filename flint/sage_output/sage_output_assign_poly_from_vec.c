#include <stdio.h>
#include <flint/nmod_poly.h>

/*---------------------------------------------------------------*/
/* prints s = U(P)                                               */
/*---------------------------------------------------------------*/
void sage_output_assign_poly_from_vec(mp_srcptr ell, long m, char *s){
  printf("%s = ", s);
  if (m == 0){
    printf("U(0)\n");
    return;
  }
    
  long i;
  for (i = 0; i < m-1; i++)
    printf("%lu*x^%ld+", ell[i], i);
  printf("%lu*x^%ld", ell[m-1], m-1);
  printf("\n");
}
