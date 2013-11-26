#include <stdio.h>
#include <flint/nmod_poly.h>
#include "sage_output.h"

/*---------------------------------------------------------------*/
/* prints s = P                                                  */
/*---------------------------------------------------------------*/
void sage_output_assign_poly(nmod_poly_t P, char *s){
  printf("%s = ", s);
  sage_output_print_poly(P);
  printf("\n");
}
