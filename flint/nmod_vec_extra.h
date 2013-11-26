#ifndef NMOD_VEC_EXTRA_H
#define NMOD_VEC_EXTRA_H

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*     A few extra functionalities for vectors over Fp                    */
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* random dense vector with n terms                                      */
/*------------------------------------------------------------------------*/
void nmod_vec_rand_dense(mp_ptr v, flint_rand_t state, long n, nmod_t mod);

/*----------------------------------------------------------------------*/
/* replaces v[i] by 1/v[i] for i=0..n-1                                 */
/* two algorithms; montgomery is usually better                         */
/* the inline functions chooses between the two                         */
/*----------------------------------------------------------------------*/

void _nmod_vec_invert_montgomery(mp_ptr v, long n, nmod_t mod);

void _nmod_vec_invert_naive(mp_ptr v, long n, nmod_t mod);

static inline void _nmod_vec_invert(mp_ptr v, long n, nmod_t mod){
  if (n < 10)
    _nmod_vec_invert_naive(v, n, mod);
  else
    _nmod_vec_invert_montgomery(v, n, mod);
}

#endif
