/*
 * fq_nmod_poly_eval.cpp
 *
 *  Created on: Dec 5, 2013
 *      Author: javad
 */

#include "fq_nmod_poly_eval.h"
#include <iostream>

using namespace std;

fq_nmod_poly_t ** fq_nmodPolyEval::init_subproduct_tree(slong num_points) {

	slong k = n_flog(num_points, 2) + 1;
	slong length = num_points;

	fq_nmod_poly_t ** tree = new fq_nmod_poly_t*[k];

	for (slong i = 0; i < k; i++) {
		tree[i] = new fq_nmod_poly_t[length];

		for (slong j = 0; j < length; j++)
			fq_nmod_poly_init(tree[i][j], ctx);

		length /= 2;
	}

	return tree;
}

void fq_nmodPolyEval::clear_subproduct_tree(fq_nmod_poly_t **tree, slong num_points) {

	slong k = n_flog(num_points, 2) + 1;
	slong length = num_points;

	for (slong i = 0; i < k; i++) {
		for (slong j = 0; j < length; j++)
			fq_nmod_poly_clear(tree[i][j], ctx);

		delete[] tree[i];
		length /= 2;
	}
}

void fq_nmodPolyEval::build_subproduct_tree(fq_nmod_poly_t **tree, const fq_nmod_t *points, slong num_points) {

	slong k = n_flog(num_points, 2) + 1;
	slong length = num_points;

	fq_nmod_t one;
	fq_nmod_t temp_coeff;

	fq_nmod_init(one, ctx);
	fq_nmod_init(temp_coeff, ctx);

	fq_nmod_set_ui(one, 1, ctx);

	for (slong i = 0; i < num_points; i++) {
		fq_nmod_neg(temp_coeff, points[i], ctx);
		fq_nmod_poly_set_coeff(tree[0][i], 0, temp_coeff, ctx);
		fq_nmod_poly_set_coeff(tree[0][i], 1, one, ctx);
	}

	for (slong i = 1; i < k; i++) {
		for (slong j = 0; j < length / 2; j++) {
			fq_nmod_poly_mul(tree[i][j], tree[i - 1][2 * j], tree[i - 1][2 * j + 1], ctx);
		}

		if (length % 2 != 0)
			fq_nmod_poly_mul(tree[i][length / 2 - 1], tree[i][length / 2 - 1], tree[i - 1][length - 1], ctx);

		length /= 2;
	}

	fq_nmod_clear(one, ctx);
	fq_nmod_clear(temp_coeff, ctx);
}

void fq_nmodPolyEval::go_down_subproduct_tree(fq_nmod_t *results, const fq_nmod_poly_t f, fq_nmod_poly_t **tree, slong num_points) {

	slong k = n_flog(num_points, 2) + 1;
	slong rev_length = n_revbin(num_points, k);
	;
	slong length = 1;

	fq_nmod_poly_rem(tree[k - 1][0], f, tree[k - 1][0], ctx);
	rev_length >>= 1;

	for (slong j = k - 2; j >= 0; j--) {
		for (slong i = 0; i < length; i++) {
			fq_nmod_poly_rem(tree[j][2 * i], tree[j + 1][i], tree[j][2 * i], ctx);
			fq_nmod_poly_rem(tree[j][2 * i + 1], tree[j + 1][i], tree[j][2 * i + 1], ctx);
		}

		if (rev_length & 1)
			fq_nmod_poly_rem(tree[j][2 * length], tree[j + 1][length - 1], tree[j][2 * length], ctx);

		length = 2 * length + (rev_length & 1);
		rev_length >>= 1;
	}

	for (slong i = 0; i < num_points; i++)
		fq_nmod_poly_get_coeff(results[i], tree[0][i], 0, ctx);

}

void fq_nmodPolyEval::multipoint_eval(fq_nmod_t *results, const fq_nmod_poly_t f, const fq_nmod_t *points, slong num_points,
		const fq_nmod_ctx_t ctx) {

	this->ctx = ctx;

	fq_nmod_poly_t **tree = init_subproduct_tree(num_points);
	build_subproduct_tree(tree, points, num_points);
	go_down_subproduct_tree(results, f, tree, num_points);
	clear_subproduct_tree(tree, num_points);

	delete[] tree;
	this->ctx = NULL;
}
