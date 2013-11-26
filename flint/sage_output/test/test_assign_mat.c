#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/nmod_poly.h>

#include "util.h"
#include "sage_output.h"


/*------------------------------------------------------------*/
/* opt is not used                                            */
/*------------------------------------------------------------*/
void check(int opt){

  mp_limb_t n = 12345;
  nmod_t Zn;
  nmod_init(&Zn, n);
  sage_output_init(Zn);

  flint_rand_t state;
  flint_randinit(state);

  nmod_mat_t a;
  nmod_mat_init(a, 10, 2, n);
  nmod_mat_randfull(a, state);
  
  sage_output_assign_mat(a, "a");

  nmod_mat_clear(a);
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
