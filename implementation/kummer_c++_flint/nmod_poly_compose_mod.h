#ifndef NMOD_POLY_COMPOSE_MOD_H_
#define NMOD_POLY_COMPOSE_MOD_H_

#include <flint/nmod_poly.h>

class Nmod_poly_compose_mod {
 public:
/*------------------------------------------------------------------------*/
/* multiple modular composition                                           */
/* res[i] = polys[i](polys[len1-1]) mod poly, i=0..n-1                    */
/* polyinv = 1 /rev(poly) mod x^n (or n-1?)                               */
/* requires n < len1                                                      */
/* the entries of res need to be uninitialised                            */
/*------------------------------------------------------------------------*/
void nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(nmod_poly_struct * res,
							 const nmod_poly_struct * polys,
							 slong len1, slong n,
							 const nmod_poly_t poly,
							 const nmod_poly_t polyinv);

};


#endif
