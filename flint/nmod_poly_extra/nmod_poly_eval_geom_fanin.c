#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"


/*-----------------------------------------------*/
/* computes the polynomial                       */
/*     f = (x-1)(x-a)(x-a^2)...(x-a^(n-1))       */
/* complexity O(M(n))                            */
/*-----------------------------------------------*/
void nmod_poly_eval_geom_fanin(nmod_poly_t f, mp_limb_t a, long N){
  long j, i = 1;

  while (i <= (N>>1))
    i<<=1;
  i>>=1;
  
  nmod_poly_fit_length(f, 2);
  f->length = 2;
  nmod_t mod  = f->mod;
  mp_limb_t n = mod.n;
  f->coeffs[0] = n-1;
  f->coeffs[1] = 1;

  nmod_poly_t tmp;
  nmod_poly_init(tmp, n);

  mp_limb_t b = a, c = a, d, ib, ia;
  ib = n_invmod(b, mod.n);
  ia = ib;

  // the algorithm is a basic divide-and-conquer
  while (i >= 1){
    nmod_poly_fit_length(tmp, f->length);
    tmp->length = f->length;
    d = c;
    long k;
    for (k = 0; k < f->length; k++){
      tmp->coeffs[k] = nmod_mul(f->coeffs[k], d, mod);
      d = nmod_mul(d, ib, mod);
    }
    nmod_poly_mul(f, f, tmp);
    b = nmod_mul(b, b, mod);
    c = nmod_mul(c, c, mod);
    c = nmod_mul(c, c, mod);
    ib = nmod_mul(ib, ib, mod);
    
    // accomodate the odd cases
    j = N & i;
    if (j != 0){
      long de = f->length;
      nmod_poly_fit_length(f, de+1);
      f->length = de+1;
      f->coeffs[de] = f->coeffs[de-1];
      long k;
      for (k = de-1; k > 0; k--){
    	mp_limb_t coef = nmod_mul(f->coeffs[k], b, mod);
    	f->coeffs[k] = nmod_sub(f->coeffs[k-1], coef, mod);
      }
      mp_limb_t coef = nmod_mul(f->coeffs[0], b, mod);
      f->coeffs[0] = nmod_sub(0, coef, mod);

      c = nmod_mul(c, b, mod);
      c = nmod_mul(c, b, mod);
      c = nmod_mul(c, a, mod);
      b = nmod_mul(b, a, mod);
      ib = nmod_mul(ib, ia, mod);
    }
    i >>= 1;
  }

  nmod_poly_clear(tmp);
}
