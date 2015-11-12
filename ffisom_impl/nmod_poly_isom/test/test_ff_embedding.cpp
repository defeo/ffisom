
#include "ff_embedding.h"
#include <iostream>
#include "nmod_min_poly.h"
#include <flint/profiler.h>

using namespace std;


void test_build_embedding(slong m, slong n, slong characteristic) {

	cout << "characteristic: " << characteristic << "\n";
	cout << "extension degrees: " << m << ", " << n << "\n";

	mp_limb_t p = characteristic;

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t f1, f2, g1, g2;
	nmod_poly_init(f1, p);
	nmod_poly_init(f2, p);
	nmod_poly_init(g1, p);
	nmod_poly_init(g2, p);

	cout << "building finite fields...\n";
	nmod_poly_randtest_monic_irreducible(f1, state, m + 1);
	nmod_poly_randtest_monic_irreducible(f2, state, n + 1);

	cout << "building an embedding...\n";
	timeit_t time;
	timeit_start(time);
	FFEmbedding ffEmbedding(f1, f2);
	ffEmbedding.compute_generators(g1, g2);
	ffEmbedding.build_embedding(g1, g2);
	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	cout << "testing the embedding... ";
	nmod_poly_t temp;
	nmod_poly_init(temp, p);
	// set temp = x
	nmod_poly_set_coeff_ui(temp, 1, 1);
	// compute the image of x
	ffEmbedding.compute_image(temp, temp);

	// the image of x should be a root of f1 in F_p[X]/(f2)
	nmod_poly_compose_mod(temp, f1, temp, f2);
	if (nmod_poly_is_zero(temp))
		cout << "ok\n";
	else
		cout << "ooops\n";

	flint_randclear(state);
	nmod_poly_clear(f1);
	nmod_poly_clear(f2);
	nmod_poly_clear(g1);
	nmod_poly_clear(g2);
	nmod_poly_clear(temp);
}


int main() {
	
	flint_rand_t state;
	flint_randinit(state);

	for (slong i = 10; i < 20; i++) {
		slong prime = n_nth_prime(2 + n_randint(state, i));
		slong m = 5 + n_randint(state, i * 30);
		slong n = m * (1 + n_randint(state, 10));
		
		// clear powers of 2
		while (m % 2 == 0)
			m /= 2;
		
		// clear powers of prime
		while (m % prime == 0)
			m /= prime;

		if (m != 1) {
			test_build_embedding(m, n, prime);
			cout << "--------------------------------\n";
		}
	}
	
	flint_randclear(state);

	return 0;
}