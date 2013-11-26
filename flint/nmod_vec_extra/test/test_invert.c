#include <flint/nmod_vec.h>

#include "util.h"
#include "sage_output.h"
#include "nmod_vec_extra.h"

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

  for (i = 1; i < 1000; i+=3){
    mp_ptr v = _nmod_vec_init(i);
    mp_ptr z = _nmod_vec_init(i);
    mp_ptr w = _nmod_vec_init(i);
    mp_ptr u = _nmod_vec_init(i);

    long j;
    for (j = 0; j < i; j++){
      mp_limb_t tmp = 0;
      do
	tmp = n_randtest_not_zero(state) % n;
      while (tmp == 0);
      v[j] = tmp;
      z[j] = tmp;
      w[j] = tmp;
      u[j] = tmp;
    }

    if (opt == 1){
      _nmod_vec_invert_montgomery(w, i, mod);
      _nmod_vec_invert_naive(u, i, mod);
      _nmod_vec_invert_naive(z, i, mod);

      sage_output_init(mod);
      sage_output_assign_vec(v, i, "v");
      sage_output_assign_vec(w, i, "w");
      sage_output_assign_vec(u, i, "u");
      sage_output_assign_vec(z, i, "z");
      printf("print w == z and w == u and w == [1/x for x in v]\n");
    }
    else{
      double t1, t2, t3;
      t1 = util_gettime();
      for (j = 0; j < 10000; j++)
	_nmod_vec_invert(w, i, mod);
      t1 = util_gettime() - t1;

      t2 = util_gettime();
      for (j = 0; j < 10000; j++)
	_nmod_vec_invert_montgomery(w, i, mod);
      t2 = util_gettime() - t2;

      t3 = util_gettime();
      for (j = 0; j < 10000; j++)
	_nmod_vec_invert_naive(w, i, mod);
      t3 = util_gettime() - t3;

      printf("%lu %f %f %f\n", i, t1, t2, t3);
    }

    _nmod_vec_clear(v);
    _nmod_vec_clear(z);
    _nmod_vec_clear(w);
    _nmod_vec_clear(u);

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
