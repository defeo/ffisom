

#include "ff_isom_artin_schreier.h"
#include <flint/nmod_poly.h>

void FFIsomArtinSchreier::compute_hilbert_90_expression(nmod_poly_t xi, nmod_poly_t beta_a, 
		nmod_poly_t beta_theta, nmod_poly_t alpha, const nmod_poly_t xi_init, 
		const nmod_poly_t a, const nmod_poly_t theta, const nmod_poly_t modulus, 
		slong n){

	if (n == 1) {
		nmod_poly_set(xi, xi_init);
		nmod_poly_set(beta_theta, theta);
		nmod_poly_set(beta_a, a);
		
		nmod_poly_compose_mod(alpha, theta, xi_init, modulus);
		nmod_poly_mulmod(alpha, alpha, a, modulus);
		
		return;
	}
	
	nmod_poly_t temp_xi;
	nmod_poly_t temp_beta_a;
	nmod_poly_t temp_beta_theta;
	nmod_poly_t temp_alpha;
	nmod_poly_t temp;
	
	nmod_poly_init(temp_xi, modulus->mod.n);
	nmod_poly_init(temp_beta_a, modulus->mod.n);
	nmod_poly_init(temp_beta_theta, modulus->mod.n);
	nmod_poly_init(temp_alpha, modulus->mod.n);
	nmod_poly_init(temp, modulus->mod.n);
	
	if ((n & 1) == 0) {
		
		compute_hilbert_90_expression(temp_xi, temp_beta_a, temp_beta_theta, 
				temp_alpha, xi_init, a, theta, modulus, n / 2);
		
		nmod_poly_compose_mod(xi, temp_xi, temp_xi, modulus);
		
		nmod_poly_compose_mod(beta_a, temp_beta_a, temp_xi, modulus);
		nmod_poly_add(beta_a, beta_a, temp_beta_a);	
		
		nmod_poly_compose_mod(beta_theta, temp_beta_theta, temp_xi, modulus);
		// need this for computing alpha
		nmod_poly_set(temp, beta_theta);
		nmod_poly_add(beta_theta, beta_theta, temp_beta_theta);
		
		nmod_poly_compose_mod(temp, temp, xi_init, modulus);
		nmod_poly_mulmod(temp, temp, temp_beta_a, modulus);
		nmod_poly_compose_mod(alpha, temp_alpha, temp_xi, modulus);
		nmod_poly_add(alpha, alpha, temp_alpha);
		nmod_poly_add(alpha, alpha, temp);
		
	} else {
		
		compute_hilbert_90_expression(temp_xi, temp_beta_a, temp_beta_theta, 
				temp_alpha, xi_init, a, theta, modulus, n - 1);
		
		nmod_poly_compose_mod(xi, temp_xi, xi_init, modulus);
		
		nmod_poly_compose_mod(beta_a, temp_beta_a, xi_init, modulus);
		nmod_poly_add(beta_a, beta_a, a);
		
		nmod_poly_compose_mod(beta_theta, temp_beta_theta, xi_init, modulus);
		// need this for computing alpha
		nmod_poly_set(temp, beta_theta);
		nmod_poly_add(beta_theta, beta_theta, theta);
		
		nmod_poly_compose_mod(temp, temp, xi_init, modulus);
		nmod_poly_mulmod(temp, temp, a, modulus);
		nmod_poly_compose_mod(alpha, temp_alpha, xi_init, modulus);
		// compute alpha1, and add it
		nmod_poly_compose_mod(temp_alpha, theta, xi_init, modulus);
		nmod_poly_mulmod(temp_alpha, temp_alpha, a, modulus);
		nmod_poly_add(alpha, alpha, temp_alpha);
		nmod_poly_add(alpha, alpha, temp);
	}
	
	nmod_poly_clear(temp_xi);
	nmod_poly_clear(temp_beta_a);
	nmod_poly_clear(temp_beta_theta);
	nmod_poly_clear(temp_alpha);
	nmod_poly_clear(temp);	
}

void FFIsomArtinSchreier::compute_hilbert_90_solution(nmod_poly_t result, const nmod_poly_t a,
		const nmod_poly_t xi_init, const nmod_poly_t modulus) {
	
	slong r = nmod_poly_degree(modulus);
	nmod_poly_t xi;
	nmod_poly_t beta_a;
	nmod_poly_t beta_theta;
	nmod_poly_t alpha;
	nmod_poly_t theta;
	
	nmod_poly_init(xi, modulus->mod.n);
	nmod_poly_init(beta_a, modulus->mod.n);
	nmod_poly_init(beta_theta, modulus->mod.n);
	nmod_poly_init(alpha, modulus->mod.n);
	nmod_poly_init(theta, modulus->mod.n);
	
	flint_rand_t state;
	flint_randinit(state);
	
	while (true) {
		nmod_poly_randtest_not_zero(theta, state, r);
		compute_hilbert_90_expression(xi, beta_a, beta_theta, alpha, xi_init, 
				a, theta, modulus, r - 1);
		
		// compute the trace of theta
		nmod_poly_compose_mod(theta, theta, xi, modulus);
		nmod_poly_add(beta_theta, beta_theta, theta);
		
		if (!nmod_poly_is_zero(beta_theta))
			break;
	}
	
	ulong theta_inverse = nmod_poly_get_coeff_ui(beta_theta, 0);
	theta_inverse = n_invmod(theta_inverse, modulus->mod.n);
	nmod_poly_scalar_mul_nmod(alpha, alpha, theta_inverse);
	nmod_poly_set(result, alpha);
	
	nmod_poly_clear(xi);
	nmod_poly_clear(beta_a);
	nmod_poly_clear(beta_theta);
	nmod_poly_clear(alpha);	
	nmod_poly_clear(theta);
	flint_randclear(state);
}


void FFIsomArtinSchreier::compute_generator(nmod_poly_t result, const nmod_poly_t modulus,
		slong exponent){
	
	nmod_poly_t const_coeff;
	nmod_poly_init(const_coeff, modulus->mod.n);
	nmod_poly_one(const_coeff);
	
	nmod_poly_t root;
	nmod_poly_init(root, modulus->mod.n);
	nmod_poly_one(root);
	
	nmod_poly_t temp;
	nmod_poly_init(temp, modulus->mod.n);
	
	nmod_poly_t xi_init;
	nmod_poly_init(xi_init, modulus->mod.n);
	nmod_poly_set_coeff_ui(xi_init, 1, 1);
	nmod_poly_powmod_ui_binexp(xi_init, xi_init, modulus->mod.n, modulus);
	
	for (slong i = 0; i < exponent; i++) {
		nmod_poly_add(temp, root, const_coeff);
		nmod_poly_mulmod(temp, temp, const_coeff, modulus);
		nmod_poly_invmod(root, root, modulus);
		nmod_poly_mulmod(const_coeff, temp, root, modulus);
		
		compute_hilbert_90_solution(root, const_coeff, xi_init, modulus);
	}

	nmod_poly_set(result, root);
	
	nmod_poly_clear(const_coeff);
	nmod_poly_clear(root);
	nmod_poly_clear(temp);
	nmod_poly_clear(xi_init);
}

void FFIsomArtinSchreier::compute_generators(nmod_poly_t g1, nmod_poly_t g2){
	slong r = nmod_poly_degree(modulus1);
	slong exponent = 0;
	while (r != 1) {
		exponent += 1;
		r /= modulus1->mod.n;
	}
	
	compute_generator(g1, modulus1, exponent);
	compute_generator(g2, modulus2, exponent);
}

FFIsomArtinSchreier::FFIsomArtinSchreier(const nmod_poly_t f1, const nmod_poly_t f2){
	nmod_poly_init(modulus1, f1->mod.n);
	nmod_poly_init(modulus2, f2->mod.n);
	nmod_poly_set(modulus1, f1);
	nmod_poly_set(modulus2, f2);
}

FFIsomArtinSchreier::~FFIsomArtinSchreier(){
	nmod_poly_clear(modulus1);
	nmod_poly_clear(modulus2);
}