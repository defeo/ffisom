#ifndef NMOD_POLY_AUTOMORPHISM_EVALUATION_H_
#define NMOD_POLY_AUTOMORPHISM_EVALUATION_H_

#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>

class Nmod_poly_automorphism_evaluation {
public:

  void _automorphism_evaluation_compose(mp_ptr res, 
					mp_srcptr a, slong len_a,
					mp_srcptr g, 
					mp_srcptr f, slong len_f,
					mp_srcptr f_inv, slong len_f_inv, nmod_t mod);

  void compose_naive(nmod_poly_t res, const nmod_poly_t A, const nmod_poly_t g, const nmod_poly_t f, const nmod_poly_t f_inv);
  void compose(nmod_poly_t res, const nmod_poly_t A, const nmod_poly_t g, const nmod_poly_t f, const nmod_poly_t f_inv);
};

#endif
