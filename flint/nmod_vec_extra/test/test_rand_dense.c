#include <flint/nmod_vec.h>

#include "util.h"
#include "sage_output.h"
#include "nmod_vec_extra.h"

/*------------------------------------------------------------*/
/* if opt = 1, runs a check                                   */
/* else, does nothing                                         */
/*------------------------------------------------------------*/
void check(int opt){
  long i;
  flint_rand_t state;
  flint_randinit(state);

  mp_limb_t n = 65537;
  nmod_t mod;
  nmod_init(&mod, n);

  for (i = 1; i < 1000; i+=3){
    mp_ptr v = _nmod_vec_init(i);
    nmod_vec_rand_dense(v, state, i, mod);

    if (opt == 1){
      sage_output_init(mod);
      sage_output_assign_vec(v, i, "v");
    }
    _nmod_vec_clear(v);
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
