#include <iostream>
#include "fq_nmod_poly_eval.h"
#include <flint/profiler.h>

using namespace std;

void test_eval(slong degree, slong num_points) {
	cout << "degree: " << degree << "\n";
	cout << "points: " << num_points << "\n";

	mp_limb_t p = 9001;

	flint_rand_t state;
	flint_randinit(state);

	nmod_poly_t f;
	nmod_poly_init(f, p);
	nmod_poly_randtest_monic_irreducible(f, state, degree);

	fq_nmod_ctx_t ctx;
	fq_nmod_ctx_init_modulus(ctx, f, "x");

	fq_nmod_poly_t g;
	fq_nmod_poly_init(g, ctx);
	fq_nmod_poly_randtest(g, state, degree, ctx);

	fq_nmod_t *points = new fq_nmod_t[num_points];
	fq_nmod_t *results = new fq_nmod_t[num_points];

	for (slong i = 0; i < num_points; i++) {
		fq_nmod_init(points[i], ctx);
		fq_nmod_init(results[i], ctx);
		fq_nmod_randtest(points[i], state, ctx);
	}

	timeit_t time;
	timeit_start(time);

	fq_nmodPolyEval fq_nmodPolyEval;
	fq_nmodPolyEval.multipoint_eval(results, g, points, num_points, ctx);

	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	fq_nmod_t temp;
	fq_nmod_init(temp, ctx);

	bool oops = false;
	for (slong i = 0; i < num_points; i++) {
		fq_nmod_poly_evaluate_fq_nmod(temp, g, points[i], ctx);
		if (!fq_nmod_equal(temp, results[i], ctx)) {
			oops = true;
			break;
		}
	}

	if (oops)
		cout << "oops\n";
	else
		cout << "ok\n";

	flint_randclear(state);
	nmod_poly_clear(f);
	fq_nmod_poly_clear(g, ctx);

	for (slong i = 0; i < num_points; i++) {
		fq_nmod_clear(points[i], ctx);
		fq_nmod_clear(results[i], ctx);
	}

	delete[] points;
	delete[] results;
	fq_nmod_clear(temp, ctx);
	fq_nmod_ctx_clear(ctx);
}

int main() {

	for (slong i = 100; i < 120; i++) {
		test_eval(i, i / 2);
		cout << "----------------------\n";
	}

	return 0;
}

