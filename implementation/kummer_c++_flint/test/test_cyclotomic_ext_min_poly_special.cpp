
#include <flint/fq_nmod_poly.h>
#include <flint/profiler.h>
#include <iostream>


#include "nmod_poly_build_irred.h"

using namespace std;


void test_cyclotomic_ext_min_poly_special(slong r, slong s, slong p) {

  cout << "r=" << r << endl;
  cout << "s=" << s << endl;
  cout << "p=" << p << endl;

  flint_rand_t state;
  flint_randinit(state);

  nmod_poly_t f;
  nmod_poly_init(f, p);
  nmod_poly_randtest_monic_irreducible(f, state, s + 1);

  fq_nmod_ctx_t ctx;
  fq_nmod_ctx_init_modulus(ctx, f, "z"); 

  fq_nmod_t a;
  fq_nmod_init(a, ctx);
  fq_nmod_randtest(a, state, ctx);

  fq_nmod_poly_t elt;
  fq_nmod_poly_init(elt, ctx);
  fq_nmod_poly_randtest_not_zero(elt, state, r, ctx);
  nmod_poly_t F;
  nmod_poly_init(F, p);

  cyclotomic_ext_min_poly_special(F, r, a, elt, r*s, ctx);

  fq_nmod_poly_t val;
  fq_nmod_poly_init(val, ctx);
  fq_nmod_poly_t mod;
  fq_nmod_poly_init(mod, ctx);
  fq_nmod_poly_set_coeff(mod, 0, a, ctx);
  fmpz_t one;
  fmpz_init(one);
  fmpz_set_ui(one, 1);
  fq_nmod_poly_set_coeff_fmpz(mod, r, one, ctx);

  nmod_poly_t tmp;
  nmod_poly_init(tmp, p);
  cout << "deg(F)=" << nmod_poly_degree(F) << endl;
  for (long i = nmod_poly_degree(F); i >= 0; i--){
    fq_nmod_poly_mul(val, val, elt, ctx);
    fq_nmod_poly_rem(val, val, mod, ctx);
    fq_nmod_poly_get_coeff(tmp, val, 0, ctx);
    nmod_poly_set_coeff_ui(tmp, 0, _nmod_add(nmod_poly_get_coeff_ui(F, i), nmod_poly_get_coeff_ui(tmp, 0), tmp->mod));
    fq_nmod_poly_set_coeff(val, 0, tmp, ctx);
  }


  cout << "F(elt)=";
  fq_nmod_poly_print_pretty(val, "x", ctx);
  cout << endl;

  nmod_poly_clear(tmp);
  fmpz_clear(one);
  fq_nmod_poly_clear(mod, ctx);
  fq_nmod_poly_clear(val, ctx);
  nmod_poly_clear(F);
  fq_nmod_poly_clear(elt, ctx);
  fq_nmod_clear(a, ctx);
  fq_nmod_ctx_clear(ctx);
  nmod_poly_clear(f);
  flint_randclear(state);

}

int main() {

  flint_rand_t state;
  flint_randinit(state);
  
  for (slong i = 10; i < 20; i++) {
    slong r = i*i+10;
    slong s = i+2;
    slong p = 65537;
    test_cyclotomic_ext_min_poly_special(r, s, p);
  }
  
  flint_randclear(state);
  
  return 0;
}
