#ifndef NMOD_POLY_COMPOSE_MOD_H_
#define NMOD_POLY_COMPOSE_MOD_H_

#include <flint/nmod_poly.h>

using namespace std;

class Nmod_poly_compose_mod {
 public:

/*------------------------------------------------------------------------*/
/* multiple modular composition                                           */
/* res[i] = polys[i](arg) mod poly, i=0..len1-1                           */
/* polyinv = 1 /rev(poly) mod x^n (or n-1?)                               */
/* sz is the number of baby steps                                         */
/* (TODO: give it a default value?)                                       */
/*------------------------------------------------------------------------*/

void nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(nmod_poly_struct * res,
							 const nmod_poly_struct * polys,
							 slong len1) const;

void nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(mp_srcptr arg, slong len_arg,
							 mp_srcptr poly, slong len_poly,
							 mp_srcptr polyinv, slong len_poly_inv,
							 nmod_t mod_in,
							 long sz);

void nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(const nmod_poly_t arg,
							 const nmod_poly_t poly,
							 const nmod_poly_t polyinv,
							 long sz);
 
 ~Nmod_poly_compose_mod();

 private:
 long n, k, m, len_f_inv;
 mp_ptr f, f_inv, giant; 
 nmod_mat_t A;  // baby steps
 nmod_t mod;
};


#endif
