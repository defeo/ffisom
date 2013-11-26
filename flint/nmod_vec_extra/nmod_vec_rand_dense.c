#include <flint/nmod_vec.h>

/*------------------------------------------------------------------------*/
/* random dense vector with n terms                                      */
/*------------------------------------------------------------------------*/

void nmod_vec_rand_dense(mp_ptr v, flint_rand_t state, long n, nmod_t mod){
  int i;
  for (i = 0; i < n; i++)
    v[i] = n_randtest(state) % mod.n;
}
