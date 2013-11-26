#include <flint/nmod_vec.h>


/*----------------------------------------------------------------------*/
/* replaces v[i] by 1/v[i] for i=0..n-1                                 */
/*----------------------------------------------------------------------*/
void _nmod_vec_invert_montgomery(mp_ptr v, long n, nmod_t mod){

  long i;
  mp_ptr v2 = _nmod_vec_init(n);
  v2[0] = v[0];
  for (i = 1; i < n; i++)
    v2[i] = nmod_mul(v[i], v2[i-1], mod);

  mp_limb_t tmp = n_invmod(v2[n-1], mod.n);
  for (i = n-1; i > 0; i--){
    mp_limb_t bak = v[i];
    v[i] = nmod_mul(v2[i-1], tmp, mod);
    tmp = nmod_mul(tmp, bak, mod);
  }

  v[0] = tmp;
  _nmod_vec_clear(v2);
}

/*----------------------------------------------------------------------*/
/* replaces v[i] by 1/v[i] for i=0..n-1                                 */
/*----------------------------------------------------------------------*/
void _nmod_vec_invert_naive(mp_ptr v, long n, nmod_t mod){
  long i;
  for (i = 0; i < n; i++)
    v[i] = n_invmod(v[i], mod.n);
}
