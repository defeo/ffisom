
#include "ff_embedding.h"
#include "nmod_min_poly.h"
#include "ff_isom_prime_power_ext.h"
#include "ff_isom_base_change.h"
#include "ff_isom_artin_schreier.h"
#include <flint/profiler.h>

#include <iostream>
using namespace std;

/**
 * Computes $\xi_{init} = x^{p^r}$
 */
void FFEmbedding::compute_xi_init(nmod_poly_t xi_init, const nmod_poly_t modulus, slong r) {
	nmod_poly_t temp_comp1;
	nmod_poly_t temp_comp2;

	nmod_poly_init(temp_comp1, modulus->mod.n);
	nmod_poly_init(temp_comp2, modulus->mod.n);

	nmod_poly_zero(temp_comp1);
	// set temp_comp1 to x
	nmod_poly_set_coeff_ui(temp_comp1, 1, 1);
	// compute x^p
	nmod_poly_powmod_ui_binexp(temp_comp1, temp_comp1, modulus->mod.n, modulus);
	nmod_poly_set(temp_comp2, temp_comp1);

	// reverse the bits of s, needed for binary-powering
	slong bit_length = n_flog(r, 2) + 1;
	slong n = n_revbin(r, bit_length);

	// compute x^{p^r} using a binary-powering scheme
	for (slong i = 1; i < bit_length; i++) {
		n >>= 1;
		nmod_poly_compose_mod(temp_comp2, temp_comp2, temp_comp2, modulus);
		if (n & 1)
			nmod_poly_compose_mod(temp_comp2, temp_comp2, temp_comp1, modulus);
	}

	nmod_poly_set(xi_init, temp_comp2);

	nmod_poly_clear(temp_comp1);
	nmod_poly_clear(temp_comp2);
}

void FFEmbedding::compute_trace(nmod_poly_t alpha, nmod_poly_t xi, const nmod_poly_t alpha_init,
		const nmod_poly_t xi_init, const nmod_poly_t modulus, slong i) {

	if (i == 1) {
		nmod_poly_set(alpha, alpha_init);
		nmod_poly_set(xi, xi_init);
		return;
	}

	nmod_poly_t temp_alpha;
	nmod_poly_t temp_xi;

	nmod_poly_init(temp_alpha, modulus->mod.n);
	nmod_poly_init(temp_xi, modulus->mod.n);

	if (i % 2 == 0) {
		compute_trace(temp_alpha, temp_xi, alpha_init, xi_init, modulus, i / 2);
		nmod_poly_compose_mod(alpha, temp_alpha, temp_xi, modulus);
		nmod_poly_add(alpha, alpha, temp_alpha);
		nmod_poly_compose_mod(xi, temp_xi, temp_xi, modulus);
	} else {
		compute_trace(temp_alpha, temp_xi, alpha_init, xi_init, modulus, i - 1);
		nmod_poly_compose_mod(alpha, temp_alpha, xi_init, modulus);
		nmod_poly_add(alpha, alpha, alpha_init);
		nmod_poly_compose_mod(xi, temp_xi, xi_init, modulus);
	}

	nmod_poly_clear(temp_alpha);
	nmod_poly_clear(temp_xi);
}

void FFEmbedding::find_subfield(nmod_poly_t subfield_modulus, nmod_poly_t embedding_image,
		const nmod_poly_t modulus, slong degree) {

	// check if the subfield is not proper
	if (degree == nmod_poly_degree(modulus)) {
		nmod_poly_set(subfield_modulus, modulus);
		// set to x
		nmod_poly_zero(embedding_image);
		nmod_poly_set_coeff_ui(embedding_image, 1, 1);

		return;
	}

	nmod_poly_t alpha;
	nmod_poly_t xi;
	nmod_poly_t alpha_init;
	nmod_poly_t xi_init;
	nmod_poly_t min_poly;

	nmod_poly_init(alpha, modulus->mod.n);
	nmod_poly_init(xi, modulus->mod.n);
	nmod_poly_init(alpha_init, modulus->mod.n);
	nmod_poly_init(xi_init, modulus->mod.n);
	nmod_poly_init(min_poly, modulus->mod.n);

	compute_xi_init(xi_init, modulus, degree);

	flint_rand_t state;
	flint_randinit(state);

	NmodMinPoly nmodMinPoly;

	// the number of terms in the trace
	slong n = nmod_poly_degree(modulus) / degree;

//	timeit_t time;
	while (true) {
		nmod_poly_randtest(alpha_init, state, nmod_poly_degree(modulus));
//		timeit_start(time);
		compute_trace(alpha, xi, alpha_init, xi_init, modulus, n);
//		timeit_stop(time);
//		cout << "trace time: " << (double) time->wall / 1000.0 << "\n";
		if (!nmod_poly_is_zero(alpha)) {
//			timeit_start(time);
			nmodMinPoly.minimal_polynomial(min_poly, alpha, modulus);
//			timeit_stop(time);
//			cout << "minpoly time: " << (double) time->wall / 1000.0 << "\n";
			if (nmod_poly_degree(min_poly) == degree)
				break;
		}
	}

	nmod_poly_set(subfield_modulus, min_poly);
	nmod_poly_set(embedding_image, alpha);

	nmod_poly_clear(alpha);
	nmod_poly_clear(xi);
	nmod_poly_clear(alpha_init);
	nmod_poly_clear(xi_init);
	nmod_poly_clear(min_poly);
	flint_randclear(state);
}


void FFEmbedding::compute_generators(nmod_poly_t g1, nmod_poly_t g2, slong r) {
	nmod_poly_t subfield_modulus1;
	nmod_poly_t subfield_modulus2;
	nmod_poly_t subfield_embd_img1;
	nmod_poly_t subfield_embd_img2;
	nmod_poly_t subfield_gen1;
	nmod_poly_t subfield_gen2;

	nmod_poly_init(subfield_modulus1, modulus1->mod.n);
	nmod_poly_init(subfield_modulus2, modulus2->mod.n);
	nmod_poly_init(subfield_embd_img1, modulus1->mod.n);
	nmod_poly_init(subfield_embd_img2, modulus2->mod.n);
	nmod_poly_init(subfield_gen1, modulus1->mod.n);
	nmod_poly_init(subfield_gen2, modulus2->mod.n);
	
//	timeit_t time;
//	timeit_start(time);
	find_subfield(subfield_modulus1, subfield_embd_img1, modulus1, r);
	find_subfield(subfield_modulus2, subfield_embd_img2, modulus2, r);
//	timeit_stop(time);
//	cout << "subfield time: " << (double) time->wall / 1000.0 << "\n";

	// check the Artin-Schreier case
	if (r % modulus1->mod.n == 0) {

		FFIsomArtinSchreier ffIsomArtinSchreier(subfield_modulus1, subfield_modulus2);
		ffIsomArtinSchreier.compute_generators(subfield_gen1, subfield_gen2);

	} else {

		FFIsomPrimePower ffIsomPrimePower(subfield_modulus1, subfield_modulus2, this->force_algo);
		ffIsomPrimePower.compute_generators(subfield_gen1, subfield_gen2);
	}

	// compute generator for subfields in the larger fields
	nmod_poly_compose_mod(g1, subfield_gen1, subfield_embd_img1, modulus1);
	nmod_poly_compose_mod(g2, subfield_gen2, subfield_embd_img2, modulus2);

	nmod_poly_clear(subfield_modulus1);
	nmod_poly_clear(subfield_modulus2);
	nmod_poly_clear(subfield_embd_img1);
	nmod_poly_clear(subfield_embd_img2);
	nmod_poly_clear(subfield_gen1);
	nmod_poly_clear(subfield_gen2);
}


void FFEmbedding::compute_generators(nmod_poly_t g1, nmod_poly_t g2) {
	slong m = nmod_poly_degree(modulus1);

	n_factor_t factors;
	n_factor_init(&factors);
	n_factor(&factors, m, 1);

	nmod_poly_t subfield_gen1;
	nmod_poly_t subfield_gen2;
	nmod_poly_init(subfield_gen1, modulus1->mod.n);
	nmod_poly_init(subfield_gen2, modulus2->mod.n);

	nmod_poly_zero(g1);
	nmod_poly_zero(g2);

	for (slong i = 0; i < factors.num; i++) {
		slong r = n_pow(factors.p[i], factors.exp[i]);
		compute_generators(subfield_gen1, subfield_gen2, r);
		nmod_poly_add(g1, g1, subfield_gen1);
		nmod_poly_add(g2, g2, subfield_gen2);
	}

	nmod_poly_clear(subfield_gen1);
	nmod_poly_clear(subfield_gen2);
}


void FFEmbedding::build_embedding(const nmod_poly_t g1, const nmod_poly_t g2) {
	nmod_poly_t x;
	nmod_poly_init(x, modulus1->mod.n);
	nmod_poly_set_coeff_ui(x, 1, 1);

	FFIsomBaseChange ffIsomBaseChange;
	ffIsomBaseChange.change_basis(x_image, g1, x, modulus1);
	nmod_poly_compose_mod(x_image, x_image, g2, modulus2);

	nmod_poly_clear(x);
}

void FFEmbedding::get_x_image(nmod_poly_t x_image) {
	nmod_poly_set(x_image, this->x_image);
}

void FFEmbedding::compute_image(nmod_poly_t image, const nmod_poly_t f) {
	nmod_poly_compose_mod(image, f, x_image, modulus2);
}

FFEmbedding::FFEmbedding(const nmod_poly_t f1, const nmod_poly_t f2, slong force_algo) {
	nmod_poly_init(modulus1, f1->mod.n);
	nmod_poly_init(modulus2, f2->mod.n);
	nmod_poly_set(modulus1, f1);
	nmod_poly_set(modulus2, f2);
	nmod_poly_init(x_image, modulus2->mod.n);
	this->force_algo = force_algo;
}

FFEmbedding::~FFEmbedding() {
	nmod_poly_clear(modulus1);
	nmod_poly_clear(modulus2);
	nmod_poly_clear(x_image);
}








