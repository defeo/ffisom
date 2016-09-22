#include "nmod_cyclotomic_poly.h"
#include "nmod_min_poly.h"
#include "nmod_fast_irred.h"
#include <math.h>
#include <iostream>

#define DEBUG 0

using namespace std;
/*
 *
 * g^(p^k) mod cyclo_r mod factor_r
 * Same routine as in construction of cyclo poly + final reduction
 * */
void NmodFastIrred::compute_qpower(nmod_poly_t result, const nmod_poly_t g, slong k, slong r, const fq_nmod_ctx_t cyclo_ctx) {
	nmod_poly_t temp;
	nmod_poly_init(temp, g->mod.n);

	// compute p^k mod r
	ulong pk = n_powmod(g->mod.n % r, k, r);

	for (slong i = 0; i < r; i++) {
		slong b = (pk * i) % r;
		nmod_poly_set_coeff_ui(temp, b, nmod_poly_get_coeff_ui(g, i));
	}

	nmod_poly_rem(result, temp, cyclo_ctx->modulus);
	nmod_poly_clear(temp);
}

/*
 * Compute alpha^{\floor{p^m/t}} modulo a factor of the r-th cyclotomic poly
 * Using Shoup 94 Lemma 14
 */
void NmodFastIrred::compute_qmtpower(fq_nmod_t alphaAk, const fq_nmod_t alpha, slong m, slong t, slong r, const fq_nmod_ctx_t cyclo_ctx) {
  // A1 = p/t
  slong A1 = cyclo_ctx->modulus->mod.n / t;

  // alpha^A1 using binexp
  fq_nmod_t alphaA1;
  fq_nmod_init(alphaA1, cyclo_ctx);
  fq_nmod_pow_ui(alphaA1, alpha, A1, cyclo_ctx);

  fq_nmod_set(alphaAk, alphaA1, cyclo_ctx);

  fq_nmod_t alphaAkqk;
  fq_nmod_init(alphaAkqk, cyclo_ctx);

  fq_nmod_t alphaAkBk;
  fq_nmod_init(alphaAkBk, cyclo_ctx);

  fq_nmod_t alphaBk2;
  fq_nmod_init(alphaBk2, cyclo_ctx);

  ulong B1 = cyclo_ctx->modulus->mod.n % t;
  ulong Bk = B1;
  ulong Bk2;

  // bitlength of m
  slong i = 0;
  for (slong j = m; j != 0; j >>= 1)
    i++;
  // double and add approach: A_{2 k + bit}
  slong k = 1;
  for (--i; i != 0; i--) {
    // double: k <- 2*k
    // (alpha^Ak)^(q^k)
    compute_qpower(alphaAkqk, alphaAk, k, r, cyclo_ctx);
    // (alpha^Ak)^Bk
    fq_nmod_pow_ui(alphaAkBk, alphaAk, Bk, cyclo_ctx);
    // alpha^(Bk²/t)
    Bk2 = Bk*Bk; //let's assume this fits in
    fq_nmod_pow_ui(alphaBk2, alpha, Bk2/t, cyclo_ctx);
    // product
    fq_nmod_mul(alphaAk, alphaAkqk, alphaAkBk, cyclo_ctx);
    fq_nmod_mul(alphaAk, alphaAk, alphaBk2, cyclo_ctx);
    Bk = Bk2 % t;
    k = m >> (i-1);
    // add: k <- k + bit
    if (k & 1) {
      // (alpha^Ak)^(q)
      compute_qpower(alphaAkqk, alphaAk, 1, r, cyclo_ctx);
      // (alpha^A1)^Bk
      fq_nmod_pow_ui(alphaAkBk, alphaA1, Bk, cyclo_ctx);
      // alpha^(Bk*B1/t)
      Bk2 = Bk*B1;
      fq_nmod_pow_ui(alphaBk2, alpha, Bk2/t, cyclo_ctx);
      // product
      fq_nmod_mul(alphaAk, alphaAkqk, alphaAkBk, cyclo_ctx);
      fq_nmod_mul(alphaAk, alphaAk, alphaBk2, cyclo_ctx);
      Bk = Bk2 % t;
    }
  }
}

slong NmodFastIrred::test_residue(const fq_nmod_t alpha, slong r, const fq_nmod_ctx_t cyclo_ctx) {
  slong result;

  // compute alpha^(q^m/r)
  fq_nmod_t beta;
  fq_nmod_init(beta, cyclo_ctx);
#if DEBUG
  printf("alpha: ");
  fq_nmod_print_pretty(alpha, cyclo_ctx);
  printf("\n");
#endif
  compute_qmtpower(beta, alpha, fq_nmod_ctx_degree(cyclo_ctx), r, r, cyclo_ctx);
#if DEBUG
  printf("beta: ");
  fq_nmod_print_pretty(beta, cyclo_ctx);
  printf("\n");
#endif

  result = fq_nmod_is_one(beta, cyclo_ctx);

  fq_nmod_clear(beta, cyclo_ctx);

  return result;
}

/*
 * Compute psi^{q^j} in (F_p[x]/(factor_r))[y]/(y^n-\xi)
 * using precomputed xiAj = xi^{\floor{p^j/n}}
 * todo: use struct with xi, xiAi, B, Bi...
 */
void NmodFastIrred::compute_qpower_ext(fq_nmod_poly_t result, const fq_nmod_poly_t psi, slong j,
	const fq_nmod_t xiAj, const fq_nmod_t xi, slong n, slong r, const fq_nmod_ctx_t cyclo_ctx) {
	// p % n
	slong B = cyclo_ctx->modulus->mod.n % n;
	// p^j % n
	slong Bj = n_powmod(B, j, n);

        fq_nmod_t alphai;
        fq_nmod_init(alphai, cyclo_ctx);

        fq_nmod_t alphaiqj;
        fq_nmod_init(alphaiqj, cyclo_ctx);

	// to store (xi^A)^i
	fq_nmod_t xiAji;
	fq_nmod_init(xiAji, cyclo_ctx);
	fq_nmod_one(xiAji, cyclo_ctx);

	fq_nmod_t xiBin;
	fq_nmod_init(xiBin, cyclo_ctx);

	fq_nmod_t coeff;
	fq_nmod_init(coeff, cyclo_ctx);

	// B*i % n = (p^j * i) % n
	slong Bji = 0;
	{
		fq_nmod_poly_get_coeff(alphai, psi, 0, cyclo_ctx);
		compute_qpower(alphaiqj, alphai, j, r, cyclo_ctx);
		fq_nmod_poly_set_coeff(result, Bji, alphaiqj, cyclo_ctx);
	}
	for (slong i = 1; i < n; i++) {
		Bji += Bj; //overflow party?
		// update (xi^A)^i
		fq_nmod_mul(xiAji, xiAji, xiAj, cyclo_ctx);
		// get coeff
		fq_nmod_poly_get_coeff(alphai, psi, i, cyclo_ctx);
		// power up coeff
		compute_qpower(alphaiqj, alphai, j, r, cyclo_ctx);
		// compute new coeff
		fq_nmod_mul(coeff, alphaiqj, xiAji, cyclo_ctx);
		// correct coeff
		fq_nmod_pow_ui(xiBin, xi, Bji/n, cyclo_ctx);
		fq_nmod_mul(coeff, coeff, xiBin, cyclo_ctx);
		// put it in its place
		fq_nmod_poly_set_coeff(result, Bji%n, coeff, cyclo_ctx);
	}

	fq_nmod_clear(alphai, cyclo_ctx);
	fq_nmod_clear(alphaiqj, cyclo_ctx);
	fq_nmod_clear(xiAji, cyclo_ctx);
	fq_nmod_clear(xiBin, cyclo_ctx);
	fq_nmod_clear(coeff, cyclo_ctx);
}

/*
 * Recursively compute the trace from (F_p[x]/(factor_r))[y]/(y^n-\xi)
 * to F_p^n up to the i-th term: sum_{j=0}^i psi^((p^(n*j))
 * and xiAi = xi^{\floor{p^(n*i)}/n}.
 */
void NmodFastIrred::compute_trace_xi(fq_nmod_poly_t trace, const fq_nmod_poly_t psi, slong i, fq_nmod_t xiAi, const fq_nmod_t xiA1, const fq_nmod_t xi, slong n, slong r, const fq_nmod_ctx_t cyclo_ctx) {

	if (i == 1) {
		fq_nmod_poly_set(trace, psi, cyclo_ctx);
		fq_nmod_set(xiAi, xiA1, cyclo_ctx);
		return;
	}

	fq_nmod_poly_t temp1;
	fq_nmod_poly_t temp2;
	fq_nmod_poly_init(temp1, cyclo_ctx);
	fq_nmod_poly_init(temp2, cyclo_ctx);

	if (i % 2 == 0) {
		// low part of the trace
		compute_trace_xi(temp1, psi, i/2, xiAi, xiA1, xi, n, r, cyclo_ctx);
		// high part of the trace
		compute_qpower_ext(temp2, temp1, n*(i/2), xiAi, xi, n, r, cyclo_ctx);
#if DEBUG
		printf("power %ld: ", n*i/2);
		fq_nmod_poly_print_pretty(temp2, "t", cyclo_ctx);
		printf("\n");
#endif
		// sum both parts
		fq_nmod_poly_add(temp1, temp1, temp2, cyclo_ctx);
		// update xiAi: double
		// mostly like in compute_qmtpower...
		// should be cleanified and refactorisied
		slong B = cyclo_ctx->modulus->mod.n % n;
		slong B1 = n_powmod(B, n, n);
		slong Bi = n_powmod(B1, (i/2), n);
		fq_nmod_t xiA2i;
		fq_nmod_init(xiA2i, cyclo_ctx);
		compute_qpower(xiA2i, xiAi, n*(i/2), r, cyclo_ctx);
		fq_nmod_t xiAiBi;
		fq_nmod_init(xiAiBi, cyclo_ctx);
		fq_nmod_pow_ui(xiAiBi, xiAi, Bi, cyclo_ctx);
		fq_nmod_t xiB2i;
		fq_nmod_init(xiB2i, cyclo_ctx);
		fq_nmod_pow_ui(xiB2i, xi, (Bi*Bi)/n, cyclo_ctx);
		fq_nmod_mul(xiAi, xiA2i, xiAiBi, cyclo_ctx);
		fq_nmod_mul(xiAi, xiAi, xiB2i, cyclo_ctx);
		fq_nmod_clear(xiA2i, cyclo_ctx);
		fq_nmod_clear(xiAiBi, cyclo_ctx);
		fq_nmod_clear(xiB2i, cyclo_ctx);
	} else {
		compute_trace_xi(temp1, psi, i - 1, xiAi, xiA1, xi, n, r, cyclo_ctx);
		compute_qpower_ext(temp2, temp1, n*1, xiA1, xi, n, r, cyclo_ctx);
		fq_nmod_poly_add(temp1, psi, temp2, cyclo_ctx);
		// update , cyclo_ctxxiAi: plus one
		// mostly like in compute_qmtpower...
		// should be cleanified and refactorisied
		slong B = cyclo_ctx->modulus->mod.n % n;
		slong B1 = n_powmod(B, n, n);
		slong Bi = n_powmod(B1, (i-1), n);
		fq_nmod_t xiA2i;
		fq_nmod_init(xiA2i, cyclo_ctx);
		compute_qpower(xiA2i, xiAi, n*1, r, cyclo_ctx);
		fq_nmod_t xiAiBi;
		fq_nmod_init(xiAiBi, cyclo_ctx);
		fq_nmod_pow_ui(xiAiBi, xiA1, Bi, cyclo_ctx);
		fq_nmod_t xiB2i;
		fq_nmod_init(xiB2i, cyclo_ctx);
		fq_nmod_pow_ui(xiB2i, xi, (Bi*B1)/n, cyclo_ctx);
		fq_nmod_mul(xiAi, xiA2i, xiAiBi, cyclo_ctx);
		fq_nmod_mul(xiAi, xiAi, xiB2i, cyclo_ctx);
		fq_nmod_clear(xiA2i, cyclo_ctx);
		fq_nmod_clear(xiAiBi, cyclo_ctx);
		fq_nmod_clear(xiB2i, cyclo_ctx);
	}
	fq_nmod_poly_set(trace, temp1, cyclo_ctx);

#if DEBUG
	printf("trace l. %ld: ", i);
	fq_nmod_poly_print_pretty(trace, "t", cyclo_ctx);
	printf("\n");
#endif
	return;
}

/*
 * Compute the trace of x from (F_p[x]/(factor_r))[y]/(y^n-\xi) down to F_p^n.
 */
void NmodFastIrred::compute_trace(fq_nmod_poly_t trace, const fq_nmod_poly_t psi, const fq_nmod_ctx_t cyclo_ctx, slong r, slong e, const fq_nmod_t xi) {
	// n = r^e
	slong n = n_pow(r, e); // overflow party

	// xi^(p^n/n)
	fq_nmod_t xiA1;
	fq_nmod_init(xiA1, cyclo_ctx);
	compute_qmtpower(xiA1, xi, n, n, r, cyclo_ctx);
	// xi^(p^(n*i)/n)
	fq_nmod_t xiAi;
	fq_nmod_init(xiAi, cyclo_ctx);

	// compute trace recursively
	compute_trace_xi(trace, psi, fq_nmod_ctx_degree(cyclo_ctx), xiAi, xiA1, xi, n, r, cyclo_ctx);

	fq_nmod_clear(xiA1, cyclo_ctx);
	fq_nmod_clear(xiAi, cyclo_ctx);

	return;
}

/*
 * Compute irreducible polynomial of degree r^e over F_p.
 * Follows Shoup 94 Alg 13
 */
void NmodFastIrred::irred_prime_power(nmod_poly_t irred, slong r, slong e, mp_limb_t p) {
  // n = r^e
  slong n = n_pow(r, e); // overflow party

  // step 1: compute factor of cyclotomic poly
  nmod_poly_t cyclo_mod;
  nmod_poly_init(cyclo_mod, p);

  NModCyclotomicPoly nModCyclotomicPoly;
  nModCyclotomicPoly.single_irred_factor(cyclo_mod, r, p);

  fq_nmod_ctx_t cyclo_ctx;
  fq_nmod_ctx_init_modulus(cyclo_ctx, cyclo_mod, "z");

#if DEBUG
  printf("cyclo_ctx: ");
  fq_nmod_ctx_print(cyclo_ctx);
  printf("\n");
#endif

  // step 2: find r-th power non-residue
  flint_rand_t state;
  flint_randinit(state);

  fq_nmod_t xi;
  fq_nmod_init(xi, cyclo_ctx);
  fq_nmod_randtest_not_zero(xi, state, cyclo_ctx);
  while (test_residue(xi, r, cyclo_ctx)) {
    fq_nmod_randtest_not_zero(xi, state, cyclo_ctx);
  }

#if DEBUG
  printf("xi: ");
  fq_nmod_print_pretty(xi, cyclo_ctx);
  printf("\n");
#endif

  // step 3: compute trace
  fq_nmod_t one;
  fq_nmod_init(one, cyclo_ctx);
  nmod_poly_one(one);
  fq_nmod_poly_t lambda;
  fq_nmod_poly_init(lambda, cyclo_ctx);
  fq_nmod_poly_set_coeff(lambda, 1, one, cyclo_ctx);
  fq_nmod_poly_t gamma;
  fq_nmod_poly_init(gamma, cyclo_ctx);
#if DEBUG
  printf("gen: ");
  fq_nmod_poly_print_pretty(lambda, "t", cyclo_ctx);
  printf("\n");
#endif
  compute_trace(gamma, lambda, cyclo_ctx, r, e, xi);
#if DEBUG
  printf("trace: ");
  fq_nmod_poly_print_pretty(gamma, "t", cyclo_ctx);
  printf("\n");
#endif

  // step 4: compute minimal polynomial
  NmodMinPoly nmodMinPoly;
  nmodMinPoly.minimal_polynomial(irred, gamma, cyclo_ctx, n, xi); 
#if DEBUG
  printf("minpoly: ");
  nmod_poly_print_pretty(irred, "x");
  printf("\n");
#endif

  nmod_poly_clear(cyclo_mod);
  fq_nmod_ctx_clear(cyclo_ctx);
  fq_nmod_clear(xi, cyclo_ctx);
  fq_nmod_clear(one, cyclo_ctx);
  fq_nmod_poly_clear(lambda, cyclo_ctx);
  fq_nmod_poly_clear(gamma, cyclo_ctx);
}
