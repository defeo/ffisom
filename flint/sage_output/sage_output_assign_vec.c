#include <stdio.h>
#include <flint/nmod_poly.h>

/*---------------------------------------------------------------*/
/* prints v = ell (vector of length m)                           */
/*---------------------------------------------------------------*/
void sage_output_assign_vec(mp_srcptr ell, long m, char *s){
  printf("%s = [", s);
  
  if (m == 0){
    printf ("]\n");
    return;
  }

  long i;
  for (i = 0; i < m-1; i++)
    printf("k(%lu), ", ell[i]);
  printf("k(%lu)]\n", ell[m-1]);
}
