#include "nmod_cyclotomic_poly.h"
#include "nmod_min_poly.h"
#include "nmod_fast_irred.h"
#include <math.h>
#include <iostream>
#include "nmod_poly_build_irred.h"
#include "util.h"
#include <flint/profiler.h>

#define DEBUG 0

using namespace std;

void NmodFastIrred::compute_qpower(nmod_poly_t result, const nmod_poly_t g, slong k, slong r) {
	nmod_poly_t temp;
	nmod_poly_init(temp, g->mod.n);

	// compute p^k mod r
	ulong pk = n_powmod(g->mod.n % r, k, r);

	for (slong i = 0; i < r; i++) {
		slong b = (pk * i) % r;
		nmod_poly_set_coeff_ui(temp, b, nmod_poly_get_coeff_ui(g, i));
	}

	nmod_poly_set(result, temp);
	nmod_poly_clear(temp);
}

/*
 *
 * g^(p^k) mod cyclo_r mod factor_r
 * Same routine as in construction of cyclo poly + final reduction
 * */
void NmodFastIrred::compute_qpower(nmod_poly_t result, const nmod_poly_t g, slong k, slong r, const fq_nmod_ctx_t cyclo_ctx) {
	compute_qpower(result, g, k, r);

	nmod_poly_rem(result, result, cyclo_ctx->modulus);
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

/*
 * Test wheter alpha is a r-th power for r dividing q^m-1.
 */
slong NmodFastIrred::test_residue(const fq_nmod_t alpha, slong r, const fq_nmod_ctx_t cyclo_ctx) {
  slong result;

  // compute alpha^(q^m/r) where r|q^m-1
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
#if DEBUG > 1
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

#if DEBUG > 1
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
void NmodFastIrred::irred_prime_power_shoup(nmod_poly_t irred, slong r, slong e, mp_limb_t p) {
  // n = r^e
  slong n = n_pow(r, e); // overflow party

  timeit_t time;

  // step 1: compute factor of cyclotomic poly
  timeit_start(time);
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

  timeit_stop(time);
  cout << "cyclo: " << (double) time->wall / 1000.0 << "\n";

  // step 2: find r-th power non-residue
  timeit_start(time);
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

  timeit_stop(time);
  cout << "residue: " << (double) time->wall / 1000.0 << "\n";

  // step 3: compute trace
  timeit_start(time);
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

  timeit_stop(time);
  cout << "trace: " << (double) time->wall / 1000.0 << "\n";

  // step 4: compute minimal polynomial
/*
  timeit_start(time);
  NmodMinPoly nmodMinPoly;
  nmodMinPoly.minimal_polynomial(irred, gamma, cyclo_ctx, n, xi); 

#if DEBUG
  printf("minpoly: ");
  nmod_poly_print_pretty(irred, "x");
  printf("\n");
#endif

  timeit_stop(time);
  cout << "minpoly (poly): " << (double) time->wall / 1000.0 << "\n";
*/

  // step 4 bis
  fq_nmod_t mxi;
  fq_nmod_init(mxi, cyclo_ctx);
  fq_nmod_neg(mxi, xi, cyclo_ctx);

  timeit_start(time);

  cyclotomic_ext_min_poly_special(irred, n, mxi, gamma, n, cyclo_ctx);

#if DEBUG
  printf("minpoly: ");
  nmod_poly_print_pretty(irred, "x");
  printf("\n");
#endif

  timeit_stop(time);
  cout << "minpoly (matrix): " << (double) time->wall / 1000.0 << "\n";

/*
  if (!nmod_poly_equal(irred, irred2)) {
    cout << "diff minpolys\n";
    abort();
  }
*/

  fq_nmod_clear(mxi, cyclo_ctx);

  nmod_poly_clear(cyclo_mod);
  fq_nmod_ctx_clear(cyclo_ctx);
  fq_nmod_clear(xi, cyclo_ctx);
  fq_nmod_clear(one, cyclo_ctx);
  fq_nmod_poly_clear(lambda, cyclo_ctx);
  fq_nmod_poly_clear(gamma, cyclo_ctx);
}

int check_degree(ulong p, slong m, slong n, slong *o) {
  ulong q;

  if (!n_is_probabprime(m))
    return 0;

  q = p % m;
  if (q == 0)
    return 0;

  Util util;
  *o = util.compute_multiplicative_order(q, m);
  return (n_gcd((m-1)/(*o), *o) == 1);// && (n_gcd(*o/n, n) == 1);
}

void NmodFastIrred::irred_prime_power_adleman_lenstra(nmod_poly_t irred, slong r, slong e, mp_limb_t p) {
  // n = r^e
  slong n = n_pow(r, e); // overflow party

  timeit_t time;

  // step 1: find roots of unity degree
  slong o;
  slong m = n + 1;
  while (!check_degree(p, m, n, &o))
    m += n;

#if DEBUG
  printf("cyclo degree: ");
  printf("%ld", m);
  printf("\n");
#endif
  // step 2: compute factor of cyclotomic poly

  timeit_start(time);
  nmod_poly_t cyclo_mod;
  nmod_poly_init(cyclo_mod, p);
  nmod_poly_set_coeff_ui(cyclo_mod, 0, p-1);
  nmod_poly_set_coeff_ui(cyclo_mod, m, 1);
  fq_nmod_ctx_t cyclo_ctx;
  fq_nmod_ctx_init_modulus(cyclo_ctx, cyclo_mod, "z");

#if DEBUG
  printf("factor: ");
  nmod_poly_print_pretty(cyclo_mod, "x");
  printf("\n");
#endif

  timeit_stop(time);
  cout << "cyclo: " << (double) time->wall / 1000.0 << "\n";

  // step 3: find root of unity
  fq_nmod_t beta;
  fq_nmod_init(beta, cyclo_ctx);
  nmod_poly_set_coeff_ui(beta, 1, 1);

  // step 4: find multiplicative generator of Z/mZ^*
  timeit_start(time);
  slong z = 2;
  Util util;
  while ((slong) util.compute_multiplicative_order(z, m) != (m-1))
    z++;
  z = n_powmod(z, o, m);
#if DEBUG
  printf("gen zm*: ");
  printf("%ld", z);
  printf("\n");
#endif

  timeit_stop(time);
  cout << "gen zm*: " << (double) time->wall / 1000.0 << "\n";

  // step 5: compute gauss period
  timeit_start(time);
  fq_nmod_t period;
  fq_nmod_init(period, cyclo_ctx);
  fq_nmod_zero(period, cyclo_ctx);
  fq_nmod_t conj;
  fq_nmod_init(conj, cyclo_ctx);
  fq_nmod_set(conj, beta, cyclo_ctx);
  for (slong i = 0; i < (m-1)/o; i++) {
    // could write z^i = a q + b and use compute_qpower for a and BSGS for b
    fq_nmod_pow_ui(conj, conj, z, cyclo_ctx);
    fq_nmod_add(period, period, conj, cyclo_ctx);
  }

  fq_nmod_clear(beta, cyclo_ctx);

#if DEBUG
  printf("period: ");
  fq_nmod_print_pretty(period, cyclo_ctx);
  printf("\n");
#endif

  timeit_stop(time);
  cout << "period: " << (double) time->wall / 1000.0 << "\n";

  // step 6: compute trace
  timeit_start(time);
  fq_nmod_t trace;
  fq_nmod_init(trace, cyclo_ctx);
  fq_nmod_zero(trace, cyclo_ctx);
  // could be merge with step 5
  for (slong i = 0; i < o/n; i++) {
    compute_qpower(conj, period, i*n, m, cyclo_ctx);
    fq_nmod_add(trace, trace, conj, cyclo_ctx);
  }
  fq_nmod_clear(period, cyclo_ctx);

#if DEBUG
  printf("trace: ");
  fq_nmod_print_pretty(trace, cyclo_ctx);
  printf("\n");
#endif

  timeit_stop(time);
  cout << "trace: " << (double) time->wall / 1000.0 << "\n";

  // step 7: compute minimal polynomial
  timeit_start(time);
  slong prec = floor(log(n)/log(p))+1; 
#if DEBUG
  printf("prec: %ld\n",prec);
#endif
  if (prec==1) {
  slong q = p;
  nmod_t qmod = cyclo_mod->mod;
  // compute root powers
  fq_nmod_t conjpows[n];
  for (slong i = 0; i < n; i++)
    fq_nmod_init(conjpows[i], cyclo_ctx);
  fq_nmod_set(conjpows[0], trace, cyclo_ctx);
  for (slong i = 1; i < n; i++)
    fq_nmod_mul(conjpows[i], conjpows[i-1], conjpows[0], cyclo_ctx);
  // compute (minus) newton sums, stored in reverse order
  mp_limb_t *newton_sums = _nmod_vec_init(n);
  slong B = p % m;
  slong Binv = n_invmod(B, m);
  for (slong i = 0; i < n; i++) {
    mp_limb_t plus = n_mulmod2_preinv(n, nmod_poly_get_coeff_ui(conjpows[i], 0), q, qmod.ninv);
    mp_limb_t minus = 0;
    // should not depend on the initial value of k
    slong k = 1;
    for (slong j = 0; j < n; j++) {
      minus = n_addmod(minus, nmod_poly_get_coeff_ui(conjpows[i], k), q);
      k = (k*Binv) % m;
    }
    newton_sums[n-1-i] = nmod_sub(minus, plus, qmod);
  }
  for (slong i = 0; i < n; i++)
    fq_nmod_clear(conjpows[i], cyclo_ctx);
  // compute coeffs (quadratic)
  mp_limb_t *coeffs = _nmod_vec_init(n+1);
  mp_limb_t invq;
  int nlimbs;
  coeffs[0] = 1;
  coeffs[1] = newton_sums[n-1];
  for (slong i = 1; i < n; i++) {
    nlimbs = _nmod_vec_dot_bound_limbs(i+1, qmod);
    coeffs[i+1] = _nmod_vec_dot(newton_sums+(n-i-1), coeffs, i+1, qmod, nlimbs);
    invq = i+1;
    n_remove(&invq, p);
    coeffs[i+1] /= (i+1)/invq;
    invq = n_invmod(i+1, q);
    coeffs[i+1] = n_mulmod2_preinv(coeffs[i+1], invq, q, qmod.ninv);
  }
  // set coeffs
  for (slong i = n; i >= 0; i--)
    nmod_poly_set_coeff_ui(irred, i, coeffs[n-i]);
  // clean up
  _nmod_vec_clear(coeffs);
  _nmod_vec_clear(newton_sums);
  } else {
  fq_nmod_poly_t minpoly;
  fq_nmod_poly_init(minpoly, cyclo_ctx);
  fq_nmod_poly_one(minpoly, cyclo_ctx);
  fq_nmod_poly_t mincopy;
  fq_nmod_poly_init(mincopy, cyclo_ctx);
  fq_nmod_set(conj, trace, cyclo_ctx);
  for (slong i = 0; i < n; i++) {
    compute_qpower(conj, conj, 1, m);
    fq_nmod_poly_set(mincopy, minpoly, cyclo_ctx);
    fq_nmod_poly_shift_left(minpoly, mincopy, 1, cyclo_ctx);
    fq_nmod_poly_scalar_submul_fq_nmod(minpoly, mincopy, conj, cyclo_ctx);
  }
  fq_nmod_poly_clear(mincopy, cyclo_ctx);

  slong coeff;
  nmod_poly_zero(irred);
  for (slong i = 0; i <= n; i++) {
    coeff = (slong) nmod_poly_get_coeff_ui(minpoly->coeffs + i, 0) - (slong) nmod_poly_get_coeff_ui(minpoly->coeffs + i, 1);
    nmod_poly_set_coeff_ui(irred, i, coeff >= 0 ? coeff : p + coeff);
  }
  }

  fq_nmod_clear(conj, cyclo_ctx);
  fq_nmod_clear(trace, cyclo_ctx);

#if DEBUG
  printf("minpoly: ");
  nmod_poly_print_pretty(irred, "x");
  printf("\n");
#endif

  timeit_stop(time);
  cout << "minpoly: " << (double) time->wall / 1000.0 << "\n";
/*
  timeit_start(time);
  NmodMinPoly nmodminpoly;
  //nmodminpoly.minimal_polynomial(irred, trace, n, cyclo_mod);
  nmodminpoly.minimal_polynomial(irred, trace, cyclo_mod);

#if DEBUG
  printf("minpoly: ");
  nmod_poly_print_pretty(irred, "x");
  printf("\n");
#endif

  timeit_stop(time);
  cout << "minpoly: " << (double) time->wall / 1000.0 << "\n";
*/

  nmod_poly_clear(cyclo_mod);
  fq_nmod_ctx_clear(cyclo_ctx);
}

void NmodFastIrred::irred_prime_power_adleman_lenstra_factor(nmod_poly_t irred, slong r, slong e, mp_limb_t p) {
  // n = r^e
  slong n = n_pow(r, e); // overflow party

  timeit_t time;

  // step 1: find roots of unity degree
  slong o;
  slong m = n + 1;
  while (!check_degree(p, m, n, &o))
    m += n;

#if DEBUG
  printf("cyclo degree: ");
  printf("%ld", m);
  printf("\n");
#endif

  // step 2: compute factor of cyclotomic poly

  timeit_start(time);
  nmod_poly_t cyclo_mod;
  nmod_poly_init(cyclo_mod, p);
  NModCyclotomicPoly nModCyclotomicPoly;
  nModCyclotomicPoly.single_irred_factor(cyclo_mod, m, p);
  fq_nmod_ctx_t cyclo_ctx;
  fq_nmod_ctx_init_modulus(cyclo_ctx, cyclo_mod, "z");

#if DEBUG
  printf("factor: ");
  nmod_poly_print_pretty(cyclo_mod, "x");
  printf("\n");
#endif

  timeit_stop(time);
  cout << "cyclo: " << (double) time->wall / 1000.0 << "\n";

  // step 3: find root of unity
  fq_nmod_t beta;
  fq_nmod_init(beta, cyclo_ctx);
  nmod_poly_set_coeff_ui(beta, 1, 1);

  // step 6: compute trace
  timeit_start(time);
  fq_nmod_t trace;
  fq_nmod_init(trace, cyclo_ctx);
  fq_nmod_zero(trace, cyclo_ctx);
  fq_nmod_t conj;
  fq_nmod_init(conj, cyclo_ctx);
  // could be merge with step 5
  for (slong i = 0; i < o/n; i++) {
    compute_qpower(conj, beta, i*n, m, cyclo_ctx);
    fq_nmod_add(trace, trace, conj, cyclo_ctx);
  }
  fq_nmod_clear(beta, cyclo_ctx);
  fq_nmod_clear(conj, cyclo_ctx);

#if DEBUG
  printf("trace: ");
  fq_nmod_print_pretty(trace, cyclo_ctx);
  printf("\n");
#endif

  timeit_stop(time);
  cout << "trace: " << (double) time->wall / 1000.0 << "\n";

  // step 7: compute minimal polynomial
  timeit_start(time);
  NmodMinPoly nmodminpoly;
  //nmodminpoly.minimal_polynomial(irred, trace, n, cyclo_mod);
  nmodminpoly.minimal_polynomial(irred, trace, cyclo_mod);
  nmod_poly_clear(trace);

#if DEBUG
  printf("minpoly: ");
  nmod_poly_print_pretty(irred, "x");
  printf("\n");
#endif

  timeit_stop(time);
  cout << "minpoly: " << (double) time->wall / 1000.0 << "\n";

  nmod_poly_clear(cyclo_mod);
  fq_nmod_ctx_clear(cyclo_ctx);
}
