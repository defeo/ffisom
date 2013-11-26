#include <flint/nmod_vec.h>

#include "util.h"
#include "sage_output.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, runs timings                                         */
/*------------------------------------------------------------*/
void check(int opt){
  long i;
  flint_rand_t state;
  flint_randinit(state);

  mp_limb_t n = 65537;
  nmod_t mod;
  nmod_init(&mod, n);

  for (i = 1; i < 2000; i+=30){

    if (opt == 1){
      mp_limb_t root = nmod_poly_find_root(i, mod);
      sage_output_init(mod);
      printf("root = k(%ld)\n", root);
      printf("root.multiplicative_order() > %ld\n", i);
    }
    else{
      double t;
      long j;
      t = util_gettime();
      for (j = 0; j < 10000; j++)
	 nmod_poly_find_root(i, mod);
      t = util_gettime()-t;

      printf("%ld %f\n", i, t);
    }
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
