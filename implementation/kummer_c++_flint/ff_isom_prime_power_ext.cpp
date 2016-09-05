/*
 * ffIsomorphism.cpp
 *
 *  Created on: Dec 3, 2013
 *      Author: javad
 */

#include "nmod_poly_compose_mod.h"
#include "nmod_poly_automorphism_evaluation.h"
#include "ff_isom_prime_power_ext.h"
#include "fq_nmod_poly_eval.h"
#include "cyclotomic_ext_rth_root.h"
#include "ff_isom_base_change.h"
#include "nmod_cyclotomic_poly.h"
#include "util.h"
#include <iostream>
#include <flint/profiler.h>

using namespace std;

/**
 * Solve HT90 using linear algebra over F_p.
 *
 * This is not restricted to prime power degree extension.
 */
void FFIsomPrimePower::compute_semi_trace_trivial_linalg(fq_nmod_t theta, const fq_nmod_ctx_t ctx, mp_limb_t z) {
        slong r = fq_nmod_ctx_degree(ctx);
        nmod_mat_t frob_auto;
        nmod_mat_init(frob_auto, r, r, ctx->modulus->mod.n);

        fq_nmod_t temp;
        fq_nmod_init(temp, ctx);
        fq_nmod_set(temp, xi_init, ctx);

        /*
 	 * Frobenius automorphism matrix.
 	 * Compute iteratively 1, x^q, x^(2q), ...
	 */
        nmod_mat_entry(frob_auto, 0, 0) = 1;
        for (slong i = 1; i < r; i++) {
                for (slong j = 0; j < r; j++) {
                        nmod_mat_entry(frob_auto, j, i) = nmod_poly_get_coeff_ui(temp, j);
                }

                fq_nmod_mul(temp, temp, xi_init, ctx);
        }

	/*
         * Frob - z Id.
         */
        for (slong i = 0; i < r; i++)
                nmod_mat_entry(frob_auto, i, i) = nmod_sub(nmod_mat_entry(frob_auto, i, i), z, ctx->modulus->mod);

        /*
         * Kernel.
         */  
        nmod_mat_nullspace(frob_auto, frob_auto);

        /*
         * Extract theta.
         */
        for (slong i = 0; i < r; i++)
                nmod_poly_set_coeff_ui(theta, i, nmod_mat_entry(frob_auto, i, 0));

        nmod_mat_clear(frob_auto);
        fq_nmod_clear(temp, ctx);

        return;
}

/**
 * Computes the value $\delta_r = a + z^{r - 1}\sigma(a) + z^{r - 2}\sigma^2(a) + \cdots + 
 * z^{1}\sigma^{r - 1}(a)$ where $a \in \mathbb{F}_p[x][z]$, and $r$ is the degree of
 * the extension {@code ctx}.
 *
 * This is done recursively using modular composition.
 */
void FFIsomPrimePower::compute_semi_trace_trivial_modcomp(fq_nmod_t theta, const fq_nmod_t a, const fq_nmod_ctx_t ctx, mp_limb_t z) {
	fq_nmod_t xi;
	fq_nmod_init(xi, ctx);
	fq_nmod_set(delta_init_trivial, a, ctx);

	_compute_semi_trace_trivial_modcomp(theta, xi, fq_nmod_ctx_degree(ctx), ctx, z);

	fq_nmod_clear(xi, ctx);
}

/**
 * Computes the value $\delta_n = a + z^{r - 1}\sigma(a) + z^{r - 2}\sigma^2(a) + \cdots + 
 * z^{r - n + 1}\sigma^{n - 1}(a)$ where $a \in \mathbb{F}_p[x][z]$ is such that $\delta_init = a$.
 * After recursion level i, $\delta_n = a + z^{r - 1}\sigma(a) + z^{r - 2}\sigma^2(a) + \cdots + 
 * z^{r - i + 1}\sigma^{i - 1}(a)$ and $\xi_i = x^{p^i}$.
 */
void FFIsomPrimePower::_compute_semi_trace_trivial_modcomp(fq_nmod_t delta, fq_nmod_t xi, slong n, const fq_nmod_ctx_t ctx, const mp_limb_t z) {

	if (n == 1) {
		fq_nmod_set(delta, delta_init_trivial, ctx);
		fq_nmod_set(xi, xi_init, ctx);

		return;
	}

	slong z_degree = 0;
	fq_nmod_t temp_delta;
	fq_nmod_t temp_xi;
	fq_nmod_init(temp_delta, ctx);
	fq_nmod_init(temp_xi, ctx);

	if (n % 2 == 0) {

		_compute_semi_trace_trivial_modcomp(delta, xi, n / 1, ctx, z);
		fq_nmod_set(temp_delta, delta, ctx);
		fq_nmod_set(temp_xi, xi, ctx);
		z_degree = fq_nmod_ctx_degree(ctx) - n / 1;

	} else {

		_compute_semi_trace_trivial_modcomp(delta, xi, n - 1, ctx, z);
		fq_nmod_set(temp_delta, delta_init_trivial, ctx);
		fq_nmod_set(temp_xi, xi_init, ctx);
		z_degree = fq_nmod_ctx_degree(ctx) - 1;

	}

	// compute delta
	nmod_poly_compose_mod(delta, delta, temp_xi, ctx->modulus);
	fq_nmod_mul_ui(delta, delta, nmod_pow_ui(z, z_degree, ctx->modulus->mod), ctx);
	fq_nmod_add(delta, delta, temp_delta, ctx);

	//compute xi
	nmod_poly_compose_mod(xi, xi, temp_xi, ctx->modulus);

	fq_nmod_clear(temp_delta, ctx);
	fq_nmod_clear(temp_xi, ctx);
}

/**
 * Solve HT90 using linear algebra over F_q.
 *
 * This is not restricted to prime power degree extension.
 */
void FFIsomPrimePower::compute_semi_trace_linalg_cyclo(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx) {
        slong r = fq_nmod_ctx_degree(ctx);
        fq_nmod_mat_t frob_auto;
        fq_nmod_mat_init(frob_auto, r, r, cyclo_ctx);

	fq_nmod_t temp;
        fq_nmod_init(temp, ctx);
        fq_nmod_set(temp, xi_init, ctx);

	fq_nmod_t cyclo_temp;
	fq_nmod_init(cyclo_temp, cyclo_ctx);

        /*
 	 * (Frob - z Id) matrix.
 	 * Compute iteratively 1, x^q, x^(2q), ...
	 */
	nmod_poly_set_coeff_ui(cyclo_temp, 0, 1);
        fq_nmod_mat_entry_set(frob_auto, 0, 0, cyclo_temp, cyclo_ctx);
        for (slong i = 1; i < r; i++) {
                for (slong j = 0; j < r; j++) {
			nmod_poly_set_coeff_ui(cyclo_temp, 0, nmod_poly_get_coeff_ui(temp, j));
			if (j == i)
				nmod_poly_set_coeff_ui(cyclo_temp, 1, nmod_neg(-1, cyclo_ctx->mod));
                        fq_nmod_mat_entry_set(frob_auto, j, i, cyclo_temp, cyclo_ctx);
                }

                fq_nmod_mul(temp, temp, xi_init, ctx);
        }

        /*
         * Kernel.
         */
        fq_nmod_mat_nullspace(frob_auto, frob_auto, cyclo_ctx);

        /*
         * Extract theta.
         */
        fq_nmod_poly_zero(theta, ctx);
        for (slong i = 0; i < fq_nmod_ctx_degree(cyclo_ctx); i++) {
		for (slong j = 0; j < r; j++)
			nmod_poly_set_coeff_ui(temp, j, nmod_poly_get_coeff_ui(fq_nmod_mat_entry(frob_auto, j, 0), i));
		fq_nmod_poly_set_coeff(theta, i, temp, ctx);
	}

        fq_nmod_mat_clear(frob_auto, cyclo_ctx);
        fq_nmod_clear(temp, ctx);
	fq_nmod_clear(cyclo_temp, cyclo_ctx);

        return;
}

/*
 * Lift solution from F_q to F_q[z]
 *
 * Uses modular exponentiation to compute frobenii.
 */
void FFIsomPrimePower::lift_ht90_modexp(fq_nmod_poly_t theta, const fq_nmod_t a0, const fq_nmod_ctx_t ctx) {
	slong s = fq_nmod_ctx_degree(cyclo_ctx);

	fq_nmod_t temp;
	fq_nmod_init(temp, ctx);
	fq_nmod_set(temp, a0, ctx);

	// a_0
	fq_nmod_poly_set_coeff(theta, 0, a0, ctx);

        // Suboptimal, use Luca's formula.
        // a_{s-1} = -1/b_0 frob(a_0)
        mp_limb_t inv_b0 = nmod_neg(nmod_inv(nmod_poly_get_coeff_ui(cyclo_mod, 0), ctx->modulus->mod), ctx->modulus->mod);
	fq_nmod_pow_ui(temp, temp, ctx->modulus->mod.n, ctx);
	fq_nmod_mul_ui(temp, temp, inv_b0, ctx);
	fq_nmod_poly_set_coeff(theta, s-1, temp, ctx);

	if (s == 2) {
		fq_nmod_clear(temp, ctx);
		return;
	}

	// s > 2
	fq_nmod_t last;
	fq_nmod_init(last, ctx);
	fq_nmod_set(last, temp, ctx);
	fq_nmod_t add;
	fq_nmod_init(add, ctx);
	// a_i = frob(a_{i+1}) + b_i a_{s-1}
        for (slong i = s-2; i > 0; i--) {
		fq_nmod_pow_ui(temp, temp, ctx->modulus->mod.n, ctx);
		fq_nmod_mul_ui(add, last, nmod_poly_get_coeff_ui(cyclo_mod, i+1), ctx);
		fq_nmod_add(temp, temp, add, ctx);
		fq_nmod_poly_set_coeff(theta, i, temp, ctx);
        }

	fq_nmod_clear(temp, ctx);
	fq_nmod_clear(last, ctx);
	fq_nmod_clear(add, ctx);
}

/*
 * Lift solution from F_q to F_q[z]
 *
 * Uses linear algebra to compute frobenii.
 */
void FFIsomPrimePower::lift_ht90_linalg(fq_nmod_poly_t theta, const fq_nmod_t a0, const nmod_mat_t frob_auto, const fq_nmod_ctx_t ctx) {
	slong r = fq_nmod_ctx_degree(ctx);
	slong s = fq_nmod_ctx_degree(cyclo_ctx);

        nmod_mat_t a[s];
        for (slong i = 0; i < s; i++)
            nmod_mat_init(a[i], r, 1, ctx->modulus->mod.n);

        // a_0
	for (slong i = 0; i < r; i++)
		nmod_mat_entry(a[0], i, 0) = nmod_poly_get_coeff_ui(a0, i);

        // Suboptimal, use Luca's formula.
        // a_{s-1}
        mp_limb_t inv_b0 = nmod_neg(nmod_inv(nmod_poly_get_coeff_ui(cyclo_mod, 0), ctx->modulus->mod), ctx->modulus->mod);
        nmod_mat_mul(a[s-1], frob_auto, a[0]);
        nmod_mat_scalar_mul(a[s-1], a[s-1], inv_b0);

	// a_i
        for (slong i = s-2; i > 0; i--) {
                nmod_mat_mul(a[i], frob_auto, a[i+1]);
                nmod_mat_scalar_mul_add(a[i], a[i], nmod_poly_get_coeff_ui(cyclo_mod, i+1), a[s-1]);
        }

	fq_nmod_t temp;
	fq_nmod_init(temp, ctx);
	for (slong i = 0; i < s; i++) {
                for (slong j = 0; j < r; j++)
                        nmod_poly_set_coeff_ui(temp, j, nmod_mat_entry(a[i], j, 0));
		fq_nmod_poly_set_coeff(theta, i, temp, ctx);
        }
	fq_nmod_clear(temp, ctx);

        for (slong i = 0; i < s; i++)
                nmod_mat_clear(a[i]);

	return;
}

/**
 * Solve HT90 using linear algebra over F_p and lifting from F_q to F_q[z].
 */
void FFIsomPrimePower::compute_semi_trace_linalg(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx) {

	slong r = fq_nmod_ctx_degree(ctx);
	nmod_mat_t frob_auto;
	nmod_mat_init(frob_auto, r, r, ctx->modulus->mod.n);

	fq_nmod_t temp;
        fq_nmod_init(temp, ctx);
        fq_nmod_set(temp, xi_init, ctx);

        // Frobenius matrix: 1, x^p, x^2p, ..., x^(r-1)p
	nmod_mat_entry(frob_auto, 0, 0) = 1;
	for (slong i = 1; i < r; i++) {
		for (slong j = 0; j < r; j++) {
			nmod_mat_entry(frob_auto, j, i) = nmod_poly_get_coeff_ui(temp, j);
		}

		fq_nmod_mul(temp, temp, xi_init, ctx);
	}

	// Suboptimal? Custom evaluation should be s r²
        nmod_mat_t cyclo_frob;
        nmod_mat_init(cyclo_frob, r, r, ctx->modulus->mod.n);
        nmod_poly_evaluate_mat(cyclo_frob, cyclo_mod, frob_auto);

        // Kernel
	nmod_mat_nullspace(cyclo_frob, cyclo_frob);

	// a_0
	for (slong i = 0; i < r; i++)
		nmod_poly_set_coeff_ui(temp, i, nmod_mat_entry(cyclo_frob, i, 0));

	// Lift a_0 to F_q
	lift_ht90_modexp(theta, temp, ctx);

        fq_nmod_clear(temp, ctx);
	nmod_mat_clear(frob_auto);
	nmod_mat_clear(cyclo_frob);

	return;
}

void FFIsomPrimePower::compute_semi_trace_cofactor(fq_nmod_poly_t theta, const nmod_poly_t cofactor, const fq_nmod_ctx_t ctx) {
	flint_rand_t state;
	flint_randinit(state);

	fq_nmod_t a;
	fq_nmod_init(a, ctx);
	fq_nmod_t a0;
	fq_nmod_init(a0, ctx);
	fq_nmod_zero(a0, ctx);

	Nmod_poly_automorphism_evaluation eval = Nmod_poly_automorphism_evaluation();

	while (fq_nmod_is_zero(a0, ctx)) {
		fq_nmod_randtest(a, state, ctx);
		eval.compose(a0, cofactor, a, ctx->modulus, ctx->inv);
	}

	lift_ht90_modexp(theta, a0, ctx);

	fq_nmod_clear(a, ctx);
	fq_nmod_clear(a0, ctx);

	flint_randclear(state);

	return;
}

/**
 * Computes the value $\delta_r = a + z^{r - 1}\sigma(a) + z^{r - 2}\sigma^2(a) + \cdots + 
 * z^{1}\sigma^{r - 1}(a)$ where $a \in \mathbb{F}_p[x][z]$, and $r$ is the degree of
 * the extension {@code ctx}.
 *
 * This is done recursively using modular composition.
 */
void FFIsomPrimePower::compute_semi_trace_modcomp(fq_nmod_poly_t theta, const fq_nmod_poly_t a, const fq_nmod_ctx_t ctx,
		const fq_nmod_poly_t cyclo_mod_lift) {

	fq_nmod_t xi;
	fq_nmod_init(xi, ctx);

	fq_nmod_poly_set(delta_init, a, ctx);

	long delta_degree = fq_nmod_poly_degree(delta_init, ctx);
	Nmod_poly_compose_mod compose_xi_init;
	compose_xi_init.nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(xi_init, ctx->modulus, ctx->inv,  
									    n_sqrt( nmod_poly_degree(ctx->modulus) * (delta_degree+2)) + 1);

	_compute_semi_trace_modcomp(theta, xi, fq_nmod_ctx_degree(ctx), compose_xi_init, ctx, cyclo_mod_lift);

	fq_nmod_clear(xi, ctx);
}


/**
 * Computes the value $\delta_n = a + z^{r - 1}\sigma(a) + z^{r - 2}\sigma^2(a) + \cdots + 
 * z^{r - n + 1}\sigma^{n - 1}(a)$ where $a \in \mathbb{F}_p[x][z]$ is such that $\delta_init = a$.
 * After recursion level i, $\delta_n = a + z^{r - 1}\sigma(a) + z^{r - 2}\sigma^2(a) + \cdots + 
 * z^{r - i + 1}\sigma^{i - 1}(a)$ and $\xi_i = x^{p^i}$.
 */
void FFIsomPrimePower::_compute_semi_trace_modcomp(fq_nmod_poly_t delta, fq_nmod_t xi, slong n, 
						   const Nmod_poly_compose_mod & compose_xi_init, 
						   const fq_nmod_ctx_t ctx,
						   const fq_nmod_poly_t cyclo_mod_lift) {

	if (n == 1) {
		fq_nmod_poly_set(delta, delta_init, ctx);
		fq_nmod_set(xi, xi_init, ctx);
		return;
	}

	slong z_degree = 0;
	fq_nmod_poly_t temp_delta;
	fq_nmod_t temp_xi;

	fq_nmod_poly_init(temp_delta, ctx);
	fq_nmod_init(temp_xi, ctx);


	if (n % 2 == 0) {

	  _compute_semi_trace_modcomp(delta, xi, n / 2, compose_xi_init, ctx, cyclo_mod_lift);
	  fq_nmod_poly_set(temp_delta, delta, ctx);
	  fq_nmod_set(temp_xi, xi, ctx);
	  z_degree = fq_nmod_ctx_degree(ctx) - n / 2;

	  if (true){
	    long delta_degree = fq_nmod_poly_degree(delta, ctx);
	    Nmod_poly_compose_mod compose;
	    compose.nmod_poly_compose_mod_brent_kung_vec_preinv_prepare(xi, ctx->modulus, ctx->inv,  
									n_sqrt( nmod_poly_degree(ctx->modulus) * (delta_degree+2)) + 1);
	    compute_delta_and_xi(delta, xi, temp_xi, z_degree, compose, ctx, cyclo_mod_lift);
	  } // we probably won't need this anymore
	  else {
	    compute_delta(delta, temp_xi, z_degree, ctx, cyclo_mod_lift);
	    compute_xi(xi, temp_xi, ctx);
	  }
		
	} else {
	  
	  _compute_semi_trace_modcomp(delta, xi, n - 1, compose_xi_init, ctx, cyclo_mod_lift);
	  fq_nmod_poly_set(temp_delta, delta_init, ctx);
	  fq_nmod_set(temp_xi, xi_init, ctx);
	  z_degree = fq_nmod_ctx_degree(ctx) - 1;

	  // an attempt to replace the modcomp's by a q-power. 
	  // not so clear if it's useful at all, so I'll leave it out for the moment
	  if (false){
	    nmod_poly_t tmp;
	    nmod_poly_init(tmp, ctx->mod.n);
	    slong deg_delta = fq_nmod_poly_degree(delta, ctx);
	    
	    for (long i = 0; i <= deg_delta; i++){
	      fq_nmod_poly_get_coeff(tmp, delta, i, ctx);
	      nmod_poly_powmod_ui_binexp_preinv(tmp, tmp, ctx->mod.n, ctx->modulus, ctx->inv);
	      fq_nmod_poly_set_coeff(delta, i, tmp, ctx);
	    }
	    nmod_poly_powmod_ui_binexp_preinv(xi, xi, ctx->mod.n, ctx->modulus, ctx->inv);
	    nmod_poly_clear(tmp);

	    shift_delta(delta, z_degree, ctx);
	  }
	  else{
	    compute_delta_and_xi(delta, xi, temp_xi, z_degree, compose_xi_init, ctx, cyclo_mod_lift);
	  }
	}
	fq_nmod_poly_add(delta, delta, temp_delta, ctx);
	fq_nmod_poly_clear(temp_delta, ctx);
	fq_nmod_clear(temp_xi, ctx);
}


 /**
 * Given $\delta = \sum_i c_i(x)z^i$, this method computes 
 * $z^{z_degree}\delta(\xi) = z^{z_degree}\sum_i c_i(\xi)z^i$. 
 * and {@code xi} = {@code xi}({@code old_xi}).
 */
void FFIsomPrimePower::compute_delta_and_xi(fq_nmod_poly_t delta, fq_nmod_t new_xi, const fq_nmod_t xi, slong z_degree, 
					    const Nmod_poly_compose_mod & compose, 
					    const fq_nmod_ctx_t ctx, const fq_nmod_poly_t cyclo_mod_lift) {

	slong delta_degree = fq_nmod_poly_degree(delta, ctx);
	slong nn = ctx->mod.n;

	// inputs to multi-modcomp
	nmod_poly_struct* input = (nmod_poly_struct *) flint_malloc((delta_degree+2)*sizeof(nmod_poly_struct));
	for (long i = 0; i <= delta_degree; i++){
	  nmod_poly_init(input + i, nn);
	  fq_nmod_poly_get_coeff(input + i, delta, i, ctx);
	}
	nmod_poly_init(input + (delta_degree + 1), nn);
	nmod_poly_set(input + (delta_degree + 1), new_xi);

	// outputs of multi-modcomp
	nmod_poly_struct* output = (nmod_poly_struct *) flint_malloc((delta_degree+2)*sizeof(nmod_poly_struct));
	compose.nmod_poly_compose_mod_brent_kung_vec_preinv_precomp(output, input, delta_degree+2);

	for (long i = 0; i <= delta_degree; i++){
	  fq_nmod_poly_set_coeff(delta, i, output + i, ctx);
	  nmod_poly_clear(output + i);
	}
	fq_nmod_set(new_xi, output + delta_degree + 1, ctx);
	nmod_poly_clear(output + delta_degree + 1);
	flint_free(output);

	for (long i = 0; i < delta_degree+2; i++)
	  nmod_poly_clear(input + i);
	flint_free(input);

	shift_delta(delta, z_degree, ctx);
}


/**
 * Computes {@code xi} = {@code xi}({@code old_xi}).
 */
void FFIsomPrimePower::compute_xi(fq_nmod_t xi, const fq_nmod_t old_xi, const fq_nmod_ctx_t ctx) {
	nmod_poly_compose_mod(xi, xi, old_xi, ctx->modulus);
}

/**
 * Given $\delta = \sum_i c_i(x)z^i$, this method computes 
 * $z^{z_degree}\delta(\xi) = z^{z_degree}\sum_i c_i(\xi)z^i$. 
 */
void FFIsomPrimePower::compute_delta(fq_nmod_poly_t delta, const fq_nmod_t xi, slong z_degree, const fq_nmod_ctx_t ctx,
		const fq_nmod_poly_t cyclo_mod_lift) {
	// temporary elements for composition
	fq_nmod_t temp_comp;

	// temporary element for coefficients
	fq_nmod_t temp_coeff;

	fq_nmod_init(temp_comp, ctx);
	fq_nmod_init(temp_coeff, ctx);
	fq_nmod_set(temp_comp, xi, ctx);

	slong delta_degree = fq_nmod_poly_degree(delta, ctx);

	// do modular composition for each coefficient
	for (slong i = 0; i <= delta_degree; i++) {

		fq_nmod_poly_get_coeff(temp_coeff, delta, i, ctx);
		nmod_poly_compose_mod(temp_coeff, temp_coeff, temp_comp, ctx->modulus);
		fq_nmod_poly_set_coeff(delta, i, temp_coeff, ctx);
	}

	// compute z^{z_degree}delta
	// TODO slightly better complexity can be obtained  
	fq_nmod_poly_shift_left(delta, delta, z_degree, ctx);
	fq_nmod_poly_rem(delta, delta, cyclo_mod_lift, ctx);

	fq_nmod_clear(temp_comp, ctx);
	fq_nmod_clear(temp_coeff, ctx);
}


// compute z^{z_degree}*delta
void FFIsomPrimePower::shift_delta(fq_nmod_poly_t delta, slong z_degree, const fq_nmod_ctx_t ctx) {
  nmod_poly_t z, z_pow;
  nmod_poly_init(z, ctx->mod.n);
  nmod_poly_init(z_pow, ctx->mod.n);

  long r = nmod_poly_degree(ctx->modulus);
  long s = nmod_poly_degree(cyclo_mod);

  nmod_poly_set_coeff_ui(z, 1, 1);
  nmod_poly_powmod_ui_binexp_preinv(z_pow, z, z_degree, cyclo_mod, cyclo_ctx->inv);
  
  nmod_poly_t coeff;
  nmod_poly_init2(coeff, ctx->mod.n, s);
  for (long j = 0; j < r; j++){
    for (long i = 0; i < s; i++)
      coeff->coeffs[i] = nmod_poly_get_coeff_ui(delta->coeffs +i , j);
    coeff->length = s;
    _nmod_poly_normalise(coeff);
    nmod_poly_mulmod(coeff, coeff, z_pow, cyclo_mod);
    for (long i = 0; i < s; i++)
      nmod_poly_set_coeff_ui(delta->coeffs +i , j, coeff->coeffs[i]);
  }
  nmod_poly_clear(coeff);

  nmod_poly_clear(z);
  nmod_poly_clear(z_pow);
}

/**
 * Computes the value $\theta = \alpha + z^{r - 1}\alpha^{p} + z^{r - 2}\alpha^{p^2} + \cdots + 
 * z\alpha^{p^{(r - 1)}$ where $\alpha \in \mathbb{F}_p[x]$, and $r$ is the degree of
 * the extension {@code ctx}.
 *
 * This is done using multi point evaluation/iterated frobenius.
 */
void FFIsomPrimePower::compute_semi_trace_iterfrob(fq_nmod_poly_t theta, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx,
		const fq_nmod_poly_t cyclo_mod_lift) {

	slong degree = fq_nmod_ctx_degree(ctx);

	fq_nmod_t frobenius[degree];
	for (slong i = 0; i < degree; i++)
		fq_nmod_init(frobenius[i], ctx);

	iterated_frobenius(frobenius, alpha, ctx, fq_nmod_poly_degree(cyclo_mod_lift, ctx));

	fq_nmod_poly_zero(theta, ctx);
	fq_nmod_poly_set_coeff(theta, 0, frobenius[0], ctx);
	for (slong i = 1; i < degree; i++) {
		fq_nmod_poly_set_coeff(theta, i, frobenius[degree - i], ctx);
	}

        // Surely suboptimal. cyclo_mod_lift has coeffs in Fp (and is fixed).
	fq_nmod_poly_rem(theta, theta, cyclo_mod_lift, ctx);

	for (slong i = 0; i < degree; i++)
		fq_nmod_clear(frobenius[i], ctx);
}

/**
 * Given an elemenet {@code alpha} in the field {@code ctx}, this method computes
 * $\alpha^{p^i}$ for all $i = 0, \dots, r - 1$. The algorithm used is Algorithm 3.1
 * form Von Zur Gathen and Shoup, 1992. More precisely, this method implements the
 * case {q = p, R = ctx, t = p, m = r - 1} of the algorithm in the paper.
 */
void FFIsomPrimePower::iterated_frobenius(fq_nmod_t *result, const fq_nmod_t alpha, const fq_nmod_ctx_t ctx, slong s) {
	fq_nmodPolyEval fq_nmodPolyEval;

	fq_nmod_poly_t temp;
	fq_nmod_poly_init(temp, ctx);

	fq_nmod_zero(result[0], ctx);
	// set result[0] to x
	nmod_poly_set_coeff_ui(result[0], 1, 1);
	fq_nmod_set(result[1], xi_init, ctx);

	slong l = n_clog(fq_nmod_ctx_degree(ctx) - 1, 2);
	slong base = 0;
	slong length = 0;

	for (slong i = 1; i <= l; i++) {

		base = 1 << (i - 1);

		// build the polynomial for multipoint evaluation
		convert(temp, result[base], ctx);

		// make sure we stay in the bound
		if (2 * base < fq_nmod_ctx_degree(ctx))
			length = base;
		else
			length = fq_nmod_ctx_degree(ctx) - base - 1;

		fq_nmodPolyEval.multipoint_eval(result + base + 1, temp, result + 1, length, ctx);
	}

	// check the trivial case of alpha = x
	if (!fq_nmod_equal(result[0], alpha, ctx)) {

		// build the polynomial for multipoint evaluation of alpha
		convert(temp, alpha, ctx);

		fq_nmodPolyEval.multipoint_eval(result, temp, result, fq_nmod_ctx_degree(ctx), ctx);
	}

	fq_nmod_poly_clear(temp, ctx);
}

/**
 * Computes a semi-trace. This methods decides the proper semi-trace computation approach
 * based on the degree of the auxiliary cyclotomic extension. 
 * 
 * @param theta		the resulting semi-trace
 */
void FFIsomPrimePower::compute_semi_trace(fq_nmod_t theta, const fq_nmod_ctx_t ctx, const mp_limb_t z) {
	fq_nmod_t alpha;
	fq_nmod_init(alpha, ctx);
	nmod_poly_set_coeff_ui(alpha, 1, 1);
	// set xi_init to x^p
	fq_nmod_pow_ui(xi_init, alpha, ctx->modulus->mod.n, ctx);

	// use naive linear algebra for low-degree module
	if (fq_nmod_ctx_degree(ctx) < linear_alg_threshold) {
                 compute_semi_trace_trivial_linalg(theta, ctx, z);
	         fq_nmod_clear(alpha, ctx);
                 return;
        }

	flint_rand_t state;
	flint_randinit(state);

        fq_nmod_zero(theta, ctx);
	while (fq_nmod_is_zero(theta, ctx)) {
		fq_nmod_randtest(alpha, state, ctx);
		compute_semi_trace_trivial_modcomp(theta, alpha, ctx, z);
	}

	fq_nmod_clear(alpha, ctx);
	flint_randclear(state);
}


// runs all algorithms
void FFIsomPrimePower::compute_semi_trace_all(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t cyclo_mod_lift) {

	slong s = fq_nmod_ctx_degree(cyclo_ctx);
	cout << "s=" << s << endl;
	flint_rand_t state;
	fq_nmod_t alpha;
	timeit_t time;

	// computing xi_init = x^p
	timeit_start(time);
	fq_nmod_init(alpha, ctx);
	nmod_poly_set_coeff_ui(alpha, 1, 1);
	fq_nmod_pow_ui(xi_init, alpha, ctx->modulus->mod.n, ctx);
	fq_nmod_clear(alpha, ctx);
	timeit_stop(time);
	cout << "x^p: " << (double) time->wall / 1000.0 << "\n";

	// linalg
	timeit_start(time);
	compute_semi_trace_linalg(theta, ctx);
	timeit_stop(time);
	cout << "linear algebra: " << (double) time->wall / 1000.0 << "\n";

	// mod comp (case 1)
	fq_nmod_poly_t alpha_poly;
	fq_nmod_poly_init(alpha_poly, ctx);
	fq_nmod_poly_zero(theta, ctx);
	timeit_start(time);
	flint_randinit(state);
	while (fq_nmod_poly_is_zero(theta, ctx)) {
	  fq_nmod_poly_randtest_not_zero(alpha_poly, state, s, ctx);
	  compute_semi_trace_modcomp(theta, alpha_poly, ctx, cyclo_mod_lift);
	}
	flint_randclear(state);
	timeit_stop(time);
	cout << "modcomp (case 1): " << (double) time->wall / 1000.0 << "\n";
	fq_nmod_poly_clear(alpha_poly, ctx);

	// automorphism evaluation (case 2)
	timeit_start(time);
	nmod_poly_t cofactor;
	nmod_poly_init(cofactor, ctx->modulus->mod.n);
	nmod_poly_zero(cofactor);
	nmod_poly_set_coeff_ui(cofactor, ext_deg, 1);
	nmod_poly_set_coeff_ui(cofactor, 0, nmod_neg(1, ctx->modulus->mod));
	nmod_poly_div(cofactor, cofactor, cyclo_mod);
	compute_semi_trace_cofactor(theta, cofactor, ctx);
	nmod_poly_clear(cofactor);
	timeit_stop(time);
	cout << "automorphism evaluation (case 2): " << (double) time->wall / 1000.0 << "\n";

	// iterated frobenius (case 3)
	// try alpha = x first
	timeit_start(time);
	fq_nmod_init(alpha, ctx);
	nmod_poly_set_coeff_ui(alpha, 1, 1);
	compute_semi_trace_iterfrob(theta, alpha, ctx, cyclo_mod_lift);
	flint_randinit(state);
	// if the semi trace of x is zero then we try random cases
	while (fq_nmod_poly_is_zero(theta, ctx)) {
	  fq_nmod_randtest_not_zero(alpha, state, ctx);
	  compute_semi_trace_iterfrob(theta, alpha, ctx, cyclo_mod_lift);
	}
	flint_randclear(state);
	timeit_stop(time);
	cout << "iterfrob (case 3): " << (double) time->wall / 1000.0 << "\n";

}

void FFIsomPrimePower::compute_semi_trace(fq_nmod_poly_t theta, const fq_nmod_ctx_t ctx, const fq_nmod_poly_t cyclo_mod_lift) {

	slong degree = fq_nmod_ctx_degree(ctx);
	slong s = fq_nmod_ctx_degree(cyclo_ctx);

	fq_nmod_t alpha;
	fq_nmod_init(alpha, ctx);
	nmod_poly_set_coeff_ui(alpha, 1, 1);
	// computing xi_init = x^p
	fq_nmod_pow_ui(xi_init, alpha, ctx->modulus->mod.n, ctx);

	// use naive linear algebra for low-degree moduli
	if (degree*s < linear_alg_threshold) {
                compute_semi_trace_linalg(theta, ctx);
	        fq_nmod_clear(alpha, ctx);
                return;
	}

	flint_rand_t state;
	flint_randinit(state);

        if (s < multi_point_threshold) {
                fq_nmod_poly_t alpha;
                fq_nmod_poly_init(alpha, ctx);
                fq_nmod_poly_zero(theta, ctx);

		timeit_t time;
		timeit_start(time);
		while (fq_nmod_poly_is_zero(theta, ctx)) {
			fq_nmod_poly_randtest_not_zero(alpha, state, s, ctx);
			compute_semi_trace_modcomp(theta, alpha, ctx, cyclo_mod_lift);
		}
		timeit_stop(time);
		cout << "semi_trace: " << (double) time->wall / 1000.0 << "\n";

                fq_nmod_poly_clear(alpha, ctx);
	} else {
	        // try alpha = x first
		compute_semi_trace_iterfrob(theta, alpha, ctx, cyclo_mod_lift);
	        // if the semi trace of x is zero then we try random cases
		while (fq_nmod_poly_is_zero(theta, ctx)) {
			fq_nmod_randtest_not_zero(alpha, state, ctx);
			compute_semi_trace_iterfrob(theta, alpha, ctx, cyclo_mod_lift);
		}
        }

	fq_nmod_clear(alpha, ctx);
	flint_randclear(state);
}


/**
 * Coerces the element {@code value} to an element of the field {@code ctx},
 * and store it in {@code result}. It is assumed that the coefficients of {@code value}
 * are in the base field.
 */
void FFIsomPrimePower::convert(fq_nmod_t result, const fq_nmod_poly_t value, const fq_nmod_ctx_t ctx) {
	fq_nmod_t temp_coeff1;
	mp_limb_t temp_coeff2 = 0;

	fq_nmod_init(temp_coeff1, ctx);
	fq_nmod_zero(result, ctx);

	for (slong i = 0; i <= fq_nmod_poly_degree(value, ctx); i++) {
		fq_nmod_poly_get_coeff(temp_coeff1, value, i, ctx);
		temp_coeff2 = nmod_poly_get_coeff_ui(temp_coeff1, 0);
		nmod_poly_set_coeff_ui(result, i, temp_coeff2);
	}

	fq_nmod_clear(temp_coeff1, ctx);
}

/**
 * Converts the element {@code value} to a polynomial over the field {@code ctx},
 * and store it in {@code result}.
 */
void FFIsomPrimePower::convert(fq_nmod_poly_t result, const fq_nmod_t value, const fq_nmod_ctx_t ctx) {
	mp_limb_t temp_coeff = 0;

	fq_nmod_t temp;
	fq_nmod_init(temp, ctx);

	fq_nmod_poly_zero(result, ctx);

	for (slong i = 0; i <= nmod_poly_degree(value); i++) {
		temp_coeff = nmod_poly_get_coeff_ui(value, i);
		nmod_poly_set_coeff_ui(temp, 0, temp_coeff);
		fq_nmod_poly_set_coeff(result, i, temp, ctx);
	}

	fq_nmod_clear(temp, ctx);
}

/**
 * Compute the isomorphism between the two extensions of the form $\mathbb{F}_p[z][x] / (x^r - \eta)$.
 * The isomorphism is of the form $x \mapsto cx$ for some $c \in \mathbb{F}_p[z]$. 
 */
void FFIsomPrimePower::compute_middle_isomorphism(fq_nmod_t c, const fq_nmod_poly_t theta_a, const fq_nmod_poly_t theta_b,
		const fq_nmod_poly_t cyclo_mod_lift) {

	fq_nmod_poly_t eta;
	fq_nmod_poly_t cyclo_mod_lift_inv_rev;

	fq_nmod_poly_init(cyclo_mod_lift_inv_rev, ctx_1);
	fq_nmod_poly_init(eta, ctx_1);

	fq_nmod_t temp;
	fq_nmod_init(temp, cyclo_ctx);

	fq_nmod_poly_reverse(cyclo_mod_lift_inv_rev, cyclo_mod_lift, fq_nmod_poly_length(cyclo_mod_lift, ctx_1), ctx_1);
	fq_nmod_poly_inv_series_newton(cyclo_mod_lift_inv_rev, cyclo_mod_lift_inv_rev, fq_nmod_poly_length(cyclo_mod_lift, ctx_1), ctx_1);

	// compute theta_a^r
	fq_nmod_poly_powmod_ui_binexp_preinv(eta, theta_a, fq_nmod_ctx_degree(ctx_1), cyclo_mod_lift, cyclo_mod_lift_inv_rev, ctx_1);
	// now eta is an element in the cyclotomic field
	convert(c, eta, cyclo_ctx);

	// compute theta_b^r
	fq_nmod_poly_powmod_ui_binexp_preinv(eta, theta_b, fq_nmod_ctx_degree(ctx_2), cyclo_mod_lift, cyclo_mod_lift_inv_rev, ctx_2);
	convert(temp, eta, cyclo_ctx);

	// compute c = theta_a^r / theta_b^r
	fq_nmod_inv(temp, temp, cyclo_ctx);
	fq_nmod_mul(c, c, temp, cyclo_ctx);

	// compute c^{1 / r}
	CyclotomicExtRthRoot cyclotomicExtRthRoot;
	cyclotomicExtRthRoot.compute_rth_root(c, c, fq_nmod_ctx_degree(ctx_1), cyclo_ctx);

	fq_nmod_poly_clear(eta, ctx_1);
	fq_nmod_poly_clear(cyclo_mod_lift_inv_rev, ctx_1);
	fq_nmod_clear(temp, cyclo_ctx);
}

/**
 * Compute the isomorphism between the cyclotomic extensions of {@code ctx_1}, {@code ctx_2}.
 * The resulting isomorphism is $f \mapsto f_{image}$.
 */
void FFIsomPrimePower::compute_extension_isomorphism(fq_nmod_poly_t f, fq_nmod_poly_t f_image) {
	fq_nmod_poly_t cyclo_mod_lift;
	fq_nmod_poly_init(cyclo_mod_lift, ctx_1);
	convert(cyclo_mod_lift, cyclo_mod, ctx_1);

	// timeit_t time;
	// timeit_start(time);
	
	if (true){
	  compute_semi_trace_all(f, ctx_1, cyclo_mod_lift);
	  compute_semi_trace_all(f_image, ctx_2, cyclo_mod_lift);
	}
	else{
	  compute_semi_trace(f, ctx_1, cyclo_mod_lift);
	  compute_semi_trace(f_image, ctx_2, cyclo_mod_lift);
	}
	// timeit_stop(time);
	// cout << "trace time: " << (double) time->wall / 1000.0 << "\n";

	fq_nmod_t c;
	fq_nmod_poly_t c_temp;
	fq_nmod_init(c, cyclo_ctx);
	fq_nmod_poly_init(c_temp, ctx_2);

	compute_middle_isomorphism(c, f, f_image, cyclo_mod_lift);

	convert(c_temp, c, ctx_2);
	fq_nmod_poly_mulmod(f_image, f_image, c_temp, cyclo_mod_lift, ctx_2);

	fq_nmod_poly_clear(cyclo_mod_lift, ctx_1);
	fq_nmod_poly_clear(c_temp, ctx_2);
	fq_nmod_clear(c, cyclo_ctx);
	fq_nmod_ctx_clear(cyclo_ctx);
}

void FFIsomPrimePower::compute_extension_isomorphism(nmod_poly_t f, nmod_poly_t f_image) {
	nmod_poly_t temp;
	nmod_poly_init(temp, ext_char);

	compute_semi_trace(f, ctx_1, cyclo_root);
	compute_semi_trace(f_image, ctx_2, cyclo_root);

	// compute the middle isomorphism
	nmod_poly_powmod_ui_binexp(temp, f, ext_deg, ctx_1->modulus);
	mp_limb_t eta1 = nmod_poly_get_coeff_ui(temp, 0);
	nmod_poly_powmod_ui_binexp(temp, f_image, ext_deg, ctx_2->modulus);
	mp_limb_t eta2 = nmod_poly_get_coeff_ui(temp, 0);
	mp_limb_t c = nmod_div(eta1, eta2, f->mod);

	CyclotomicExtRthRoot cyclotomicExtRthRoot;
	mp_limb_t root = cyclotomicExtRthRoot.compute_rth_root(c, ext_deg, ext_char);

	nmod_poly_scalar_mul_nmod(f_image, f_image, root);

	nmod_poly_clear(temp);
}

void FFIsomPrimePower::compute_generators_trivial(nmod_poly_t g1, nmod_poly_t g2) {
	nmod_poly_t f;
	nmod_poly_t f_image;
	nmod_poly_init(f, g1->mod.n);
	nmod_poly_init(f_image, g2->mod.n);

	compute_extension_isomorphism(f, f_image);

	nmod_poly_set(g1, f);
	nmod_poly_set(g2, f_image);

	nmod_poly_clear(f);
	nmod_poly_clear(f_image);

	return;
}

void FFIsomPrimePower::compute_generators_nontriv(nmod_poly_t g1, nmod_poly_t g2) {
	fq_nmod_poly_t f;
	fq_nmod_poly_t f_image;

	fq_nmod_poly_init(f, ctx_1);
	fq_nmod_poly_init(f_image, ctx_2);

	compute_extension_isomorphism(f, f_image);

	fq_nmod_poly_get_coeff(g1, f, 0, ctx_1);
	fq_nmod_poly_get_coeff(g2, f_image, 0, ctx_2);

	fq_nmod_poly_clear(f, ctx_1);
	fq_nmod_poly_clear(f_image, ctx_2);

        return;
}

void FFIsomPrimePower::compute_generators(nmod_poly_t g1, nmod_poly_t g2) {

	if (cyclo_deg == 1)
                compute_generators_trivial(g1, g2);
        else
                compute_generators_nontriv(g1, g2);

}

/**
 * Builds a cyclotomic extension of the prime field by computing an irreducible
 * factor of the $r$-th cyclotomic polynomial.
 */
void FFIsomPrimePower::compute_cyclotomic_extension() {
	NModCyclotomicPoly nModCyclotomicPoly;
	nModCyclotomicPoly.single_irred_factor(cyclo_mod, ext_deg, ext_char);
	fq_nmod_ctx_init_modulus(cyclo_ctx, cyclo_mod, "z");
}

void FFIsomPrimePower::compute_cyclotomic_root() {
	NModCyclotomicPoly nModCyclotomicPoly;
	nModCyclotomicPoly.single_irred_factor(cyclo_mod, ext_deg, ext_char);
	cyclo_root = ext_char - nmod_poly_get_coeff_ui(cyclo_mod, 0);
}

FFIsomPrimePower::FFIsomPrimePower(const nmod_poly_t modulus1, 
				   const nmod_poly_t modulus2, slong linear_alg_threshold, slong multi_point_threshold) {
	Util util;

	this->linear_alg_threshold = linear_alg_threshold;
	this->multi_point_threshold = multi_point_threshold;

	ext_char = modulus1->mod.n;
        ext_deg = nmod_poly_degree(modulus1);

	fq_nmod_ctx_init_modulus(ctx_1, modulus1, "x");
	fq_nmod_ctx_init_modulus(ctx_2, modulus2, "x");

	fq_nmod_poly_init(delta_init, ctx_1);
	fq_nmod_init(xi_init, ctx_1);
	fq_nmod_init(delta_init_trivial, ctx_1);

	// check for the trivial cyclotomic extension case
	cyclo_deg = util.compute_multiplicative_order(ext_char, ext_deg);
	nmod_poly_init(cyclo_mod, ext_char);
        if (cyclo_deg == 1) {
		compute_cyclotomic_root();
	} else {
		compute_cyclotomic_extension();
	}

}


FFIsomPrimePower::~FFIsomPrimePower() {
	fq_nmod_poly_clear(delta_init, ctx_1);
	fq_nmod_clear(xi_init, ctx_1);
	fq_nmod_clear(delta_init_trivial, ctx_1);
	fq_nmod_ctx_clear(ctx_1);
	fq_nmod_ctx_clear(ctx_2);
	nmod_poly_clear(cyclo_mod);
	if (cyclo_deg > 1) {
	  fq_nmod_ctx_clear(cyclo_ctx);
	}
}
