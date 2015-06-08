#include <iostream>
#include "fq_poly_eval.h"
#include <flint/profiler.h>

using namespace std;

void test_eval(slong degree, slong num_points) {
	cout << "degree: " << degree << "\n";
	cout << "points: " << num_points << "\n";

	fmpz_t p;
	fmpz_init_set_ui(p, 9001);

	flint_rand_t state;
	flint_randinit(state);

	fmpz_mod_poly_t f;
	fmpz_mod_poly_init(f, p);
	fmpz_mod_poly_randtest_monic_irreducible(f, state, degree);

	fq_ctx_t ctx;
	fq_ctx_init_modulus(ctx, f, "x");

	fq_poly_t g;
	fq_poly_init(g, ctx);
	fq_poly_randtest(g, state, degree, ctx);

	fq_t *points = new fq_t[num_points];
	fq_t *results = new fq_t[num_points];

	for (slong i = 0; i < num_points; i++) {
		fq_init(points[i], ctx);
		fq_init(results[i], ctx);
		fq_randtest(points[i], state, ctx);
	}

	timeit_t time;
	timeit_start(time);
	
	FqPolyEval fqPolyEval;
	fqPolyEval.multipoint_eval(results, g, points, num_points, ctx);
	
	timeit_stop(time);
	cout << "time: " << (double) time->wall / 1000.0 << "\n";

	fq_t temp;
	fq_init(temp, ctx);

	bool oops = false;
	for (slong i = 0; i < num_points; i++) {
		fq_poly_evaluate_fq(temp, g, points[i], ctx);
		if (!fq_equal(temp, results[i], ctx)) {
			oops = true;
			break;
		}
	}

	if (oops)
		cout << "oops\n";
	else
		cout << "ok\n";

	fmpz_clear(p);
	flint_randclear(state);
	fmpz_mod_poly_clear(f);
	fq_poly_clear(g, ctx);

	for (slong i = 0; i < num_points; i++) {
		fq_clear(points[i], ctx);
		fq_clear(results[i], ctx);
	}

	delete[] points;
	delete[] results;
	fq_clear(temp, ctx);
	fq_ctx_clear(ctx);
}


int main() {
	
	for (slong i = 100; i < 120; i++){
		test_eval(i, i / 2);
		cout << "----------------------\n";
	}
	
	return 0;
}

