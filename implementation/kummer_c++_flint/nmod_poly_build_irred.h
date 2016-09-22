#ifndef NMOD_POLY_BUILD_IRRED_H_
#define NMOD_POLY_BUILD_IRRED_H_

#include <flint/nmod_poly.h>

//---------------------------------------------------
// computes F = min-poly of elt modulo <x^r+a, ctx> 
// assumes that we are in a field extension (no check is done)
// t is a bound on the degree of the minpoly we want to compute
//----------------------------------------------------
void cyclotomic_ext_min_poly_special(nmod_poly_t F, slong r, fq_nmod_t a, fq_nmod_poly_t elt, slong t, fq_nmod_ctx_t ctx);

#endif

