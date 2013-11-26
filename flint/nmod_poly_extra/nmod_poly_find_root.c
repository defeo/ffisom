#include <flint/nmod_vec.h>

/*--------------------------------------------------------*/
/* finds an element of order greater than ord modulo mod  */
/* if no such element exists, return 0                    */
/* starts the search at start                             */
/*--------------------------------------------------------*/
mp_limb_t nmod_poly_find_root_seed(long ord, long start, nmod_t mod){

  long i;
  for (i = 0; i < mod.n; i++){
    mp_limb_t loc = nmod_add(start, i, mod);
    
    if (loc == 0)
      continue;

    long done = 1;
    mp_limb_t rho = loc, tmp = 1;
    long j;
    for (j = 1; j <= ord; j++){
      tmp = nmod_mul(tmp, rho, mod); 
      if (tmp == 1)  // invariant: tmp = rho^j
	done = 0;
    }
    if (done == 1)
      return loc;
  }
  
  return 0;
}

mp_limb_t nmod_poly_find_root(long ord, nmod_t mod){
  return nmod_poly_find_root_seed(ord, 2, mod);
}
