#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/nmod_poly.h>

#include "util.h"
#include "sage_output.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* opt is not used                                            */
/*------------------------------------------------------------*/
void check(int opt){
  long i;
  flint_rand_t state;
  flint_randinit(state);

  mp_limb_t n = 65537;

  for (i = 0; i < 100; i+=10){

    nmod_poly_t P;
    nmod_poly_init(P, n);
    nmod_poly_rand_dense(P, state, i);

    sage_output_init(P->mod);
    sage_output_assign_poly(P, "P");
    nmod_poly_clear(P);

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


