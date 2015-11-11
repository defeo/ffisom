#include <iostream>
#include "ff_isom.h"
#include "util.h"
#include "nmod_min_poly.h"
#include <flint/profiler.h>

using namespace std;

void test_build_isom(slong degree, slong characteristic) {

	Util util;
	cout << "characteristic: " << characteristic << "\n";
	cout << "extension degree: " << degree << "\n";
	cout << "cyclotomic extension degree: " << util.compute_multiplicative_order(characteristic, degree) << "\n";

	mp_limb_t p = characteristic;

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t f1, f2, g1, g2;
	nmod_poly_init(f1, p);
	nmod_poly_init(f2, p);
	nmod_poly_init(g1, p);
	nmod_poly_init(g2, p);

	cout << "building finite fields...\n";
	nmod_poly_randtest_monic_irreducible(f1, state, degree + 1);
	NmodMinPoly fmpzMinPoly;
	while (nmod_poly_degree(f2) < degree) {
		nmod_poly_randtest_not_zero(f2, state, degree);
		fmpzMinPoly.minimal_polynomial(f2, f2, f1);
	}

	cout << "building an isomorphism...\n";
	timeit_t time;
	timeit_start(time);
	FFIsomorphism ffIsomorphism(f1, f2);
	ffIsomorphism.compute_generators(g1, g2);
	ffIsomorphism.build_isomorphism(g1, g2);
	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	cout << "testing the isomorphism... ";
	nmod_poly_t temp;
	nmod_poly_init(temp, p);
	// set temp = x
	nmod_poly_set_coeff_ui(temp, 1, 1);
	// compute the image of x
	ffIsomorphism.compute_isom_matrix();
	ffIsomorphism.compute_image_using_matrix(temp, temp);

	// the image of x should be a root of f1 in F_p[X]/<f2>
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

double test_build_isom1(slong degree, slong prime) {

	mp_limb_t p = prime;

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t f1, f2, g1, g2;
	nmod_poly_init(f1, p);
	nmod_poly_init(f2, p);
	nmod_poly_init(g1, p);
	nmod_poly_init(g2, p);

	nmod_poly_randtest_monic_irreducible(f1, state, degree + 1);
	NmodMinPoly fmpzMinPoly;
	while (nmod_poly_degree(f2) < degree) {
		nmod_poly_randtest_not_zero(f2, state, degree);
		fmpzMinPoly.minimal_polynomial(f2, f2, f1);
	}

	timeit_t time;
	timeit_start(time);
	FFIsomorphism ffIsomorphism(f1, f2);
	ffIsomorphism.compute_generators(g1, g2);
	ffIsomorphism.build_isomorphism(g1, g2);
	timeit_stop(time);

	flint_randclear(state);
	nmod_poly_clear(f1);
	nmod_poly_clear(f2);
	nmod_poly_clear(g1);
	nmod_poly_clear(g2);
	return (double) time->wall / 1000.0;
}

void test_prime_ext() {

	flint_rand_t state;
	flint_randinit(state);

	cout << "\n/////////////////////// random case //////////////////////////\n\n";
	//	cout << "prime  extDegree  time \n";	
	for (slong i = 30; i < 50; i++) {
		slong prime = n_nth_prime(5 + n_randint(state, i));
		slong degree = n_nth_prime(5 + n_randint(state, i));

		test_build_isom(degree, prime);
		cout << "--------------------------------\n";
		//		double time = test_build_isom1(degree, prime);
		//		cout << prime << "  " << degree << "  " << time << "\n";
	}
	flint_randclear(state);

	cout << "\n/////////////////////// no cyclotomic extension //////////////////////////\n\n";
	//	cout << "prime  extDegree  time \n";
	n_factor_t factors;
	for (slong i = 100; i < 120; i++) {//1048576
		slong prime = n_nth_prime(i);
		n_factor_init(&factors);
		n_factor(&factors, prime - 1, 1);
		slong degree = factors.p[factors.num - 1];

		test_build_isom(degree, prime);
		cout << "--------------------------------\n";
		//		double time = test_build_isom1(degree, prime);
		//		cout << prime << "  " << degree << "  " << time << "\n";
	}

	Util util;
	cout << "\n/////////////////////// large cyclotomic extension //////////////////////////\n\n";
	cout << "prime  extDegree  time \n";
	for (slong i = 20; i < 50; i++) {
		slong prime = n_nth_prime(i);
		slong degree = n_nth_prime(i - 5);
		if (util.compute_multiplicative_order(prime, degree) == n_euler_phi(degree)) {

			test_build_isom(degree, prime);
			cout << "--------------------------------\n";
			//			double time = test_build_isom1(degree, prime);
			//			cout << prime << "  " << degree << "  " << time << "\n";
		}
	}

}

void test_prime_power_ext() {

	flint_rand_t state;
	flint_randinit(state);

	slong exp = 2;
	cout << "\n/////////////////////// random case //////////////////////////\n\n";
	//	cout << "prime  extDegree  time \n";	
	for (slong i = 2; i < 10; i++) {
		slong prime = n_nth_prime(2 + n_randint(state, i));
		slong degree = n_nth_prime(2 + n_randint(state, i));
		degree = n_pow(degree, exp);

		test_build_isom(degree, prime);
		cout << "--------------------------------\n";
		//		double time = test_build_isom1(degree, prime);
		//		cout << prime << "  " << degree << "  " << time << "\n";
	}
	flint_randclear(state);

	Util util;
	cout << "\n/////////////////////// large cyclotomic extension //////////////////////////\n\n";
	cout << "prime  extDegree  time \n";
	for (slong i = 7; i < 15; i++) {
		slong prime = n_nth_prime(i);
		slong degree = n_nth_prime(i - 5);
		degree = n_pow(degree, exp);
		if (util.compute_multiplicative_order(prime, degree) == n_euler_phi(degree)) {

			test_build_isom(degree, prime);
			cout << "--------------------------------\n";
			//			double time = test_build_isom1(degree, prime);
			//			cout << prime << "  " << degree << "  " << time << "\n";
		}
	}
}

int main() {

	test_prime_power_ext();
	return 0;
}

