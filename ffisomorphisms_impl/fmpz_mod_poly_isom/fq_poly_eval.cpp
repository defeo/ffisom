/*
 * fq_poly_eval.cpp
 *
 *  Created on: Dec 5, 2013
 *      Author: javad
 */

#include "fq_poly_eval.h"
#include <iostream>

using namespace std;

fq_poly_t ** FqPolyEval::init_subproduct_tree(slong num_points) {

	slong k = n_flog(num_points, 2) + 1;
	slong length = num_points;

	fq_poly_t ** tree = new fq_poly_t*[k];

	for (slong i = 0; i < k; i++) {
		tree[i] = new fq_poly_t[length];

		for (slong j = 0; j < length; j++)
			fq_poly_init(tree[i][j], ctx);
		
		length /= 2;
	}

	return tree;
}

void FqPolyEval::clear_subproduct_tree(fq_poly_t **tree, slong num_points) {

	slong k = n_flog(num_points, 2) + 1;
	slong length = num_points;

	for (slong i = 0; i < k; i++) {
		for (slong j = 0; j < length; j++)
			fq_poly_clear(tree[i][j], ctx);

		delete[] tree[i];
		length /= 2;
	}
}

void FqPolyEval::build_subproduct_tree(fq_poly_t **tree, const fq_t *points, slong num_points) {

	slong k = n_flog(num_points, 2) + 1;
	slong length = num_points;

	fq_t one;
	fq_t temp_coeff;

	fq_init(one, ctx);
	fq_init(temp_coeff, ctx);

	fq_set_ui(one, 1, ctx);

	for (slong i = 0; i < num_points; i++) {
		fq_neg(temp_coeff, points[i], ctx);
		fq_poly_set_coeff(tree[0][i], 0, temp_coeff, ctx);
		fq_poly_set_coeff(tree[0][i], 1, one, ctx);
	}

	for (slong i = 1; i < k; i++) {
		for (slong j = 0; j < length / 2; j++) {
			fq_poly_mul(tree[i][j], tree[i - 1][2 * j], tree[i - 1][2 * j + 1], ctx);
		}
		
		if (length % 2 != 0)
			fq_poly_mul(tree[i][length / 2 - 1], tree[i][length / 2 - 1], tree[i - 1][length - 1], ctx);
		
		length /= 2;
	}

	fq_clear(one, ctx);
	fq_clear(temp_coeff, ctx);
}

void FqPolyEval::go_down_subproduct_tree(fq_t *results, const fq_poly_t f, fq_poly_t **tree, slong num_points) {

	slong k = n_flog(num_points, 2) + 1;
	slong rev_length = n_revbin(num_points, k);;
	slong length = 1;

	fq_poly_rem(tree[k - 1][0], f, tree[k - 1][0], ctx);
	rev_length >>= 1;

	for (slong j = k - 2; j >= 0; j--) {
		for (slong i = 0; i < length; i++) {
			fq_poly_rem(tree[j][2 * i], tree[j + 1][i], tree[j][2 * i], ctx);
			fq_poly_rem(tree[j][2 * i + 1], tree[j + 1][i], tree[j][2 * i + 1], ctx);
		}
		
		if (rev_length & 1)
			fq_poly_rem(tree[j][2 * length], tree[j + 1][length - 1], tree[j][2 * length], ctx);
			
		length = 2 * length + (rev_length & 1);
		rev_length >>= 1;
	}

	for (slong i = 0; i < num_points; i++)
		fq_poly_get_coeff(results[i], tree[0][i], 0, ctx);

}

void FqPolyEval::multipoint_eval(fq_t *results, const fq_poly_t f, const fq_t *points, slong num_points,
		const fq_ctx_t ctx) {

	this->ctx = ctx;

	fq_poly_t **tree = init_subproduct_tree(num_points);
	build_subproduct_tree(tree, points, num_points);
	go_down_subproduct_tree(results, f, tree, num_points);
	clear_subproduct_tree(tree, num_points);

	delete[] tree;
	this->ctx = NULL;
}
