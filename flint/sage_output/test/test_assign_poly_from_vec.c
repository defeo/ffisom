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

  mp_ptr a = _nmod_vec_init(10);
  long i;
  for (i = 0; i < 10; i++)
    a[i] = i;

  sage_output_assign_poly_from_vec(a, 10, "a");
  sage_output_assign_poly_from_vec(a, 0, "b");

  _nmod_vec_clear(a);
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
