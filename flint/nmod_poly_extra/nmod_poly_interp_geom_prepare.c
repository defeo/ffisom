#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"



/*------------------------------------------------------------------------*/
/* G = reverse(fanin)                                                     */
/* inverse_derivative[i] = 1/(fanin'(q^(2i)) q^((n-i)^2), i < n           */
/*------------------------------------------------------------------------*/
void nmod_poly_interp_geom_prepare(nmod_poly_t G, mp_ptr inverse_derivative, 
				   mp_srcptr inverse_powers_square_q, mp_limb_t q, long n){


  nmod_poly_t F;
  nmod_t mod = G->mod;
  nmod_poly_init(F, mod.n);
  mp_limb_t a = nmod_mul(q, q, mod);

  nmod_poly_eval_geom_fanin(F, a, n);
  nmod_poly_reverse(G, F, n+1);

  if (n == 1){
    inverse_derivative[0] = n_invmod(q, mod.n);
    nmod_poly_clear(F);
    return;
  }

  nmod_poly_eval_geom_derivative_fanin(inverse_derivative, a, n, mod);
  _nmod_vec_invert_montgomery(inverse_derivative, n, mod);

  // values of F' and their inverses, premultiplied by InvPowSqQ[n-i]
  long i;
  for (i = 0; i < n; i++)
    inverse_derivative[i] = nmod_mul(inverse_derivative[i], inverse_powers_square_q[n-i], mod);

  nmod_poly_clear(F);
}

