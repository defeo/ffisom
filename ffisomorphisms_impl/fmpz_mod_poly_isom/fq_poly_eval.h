/*
 * fq_poly_multippoint_eval.h
 *
 *  Created on: Dec 5, 2013
 *      Author: javad
 */

#ifndef FQ_POLY_MULTIPPOINT_EVAL_H_
#define FQ_POLY_MULTIPPOINT_EVAL_H_

#include <flint/fq_poly.h>


/**
 * This class is fast multipoint evaluation over finite fields.
 */
class FqPolyEval {

	const fq_ctx_struct *ctx;

	void build_subproduct_tree(fq_poly_t **tree, const fq_t *points, slong num_points);
	void go_down_subproduct_tree(fq_t *results, const fq_poly_t f, fq_poly_t **tree, slong num_points);
	fq_poly_t ** init_subproduct_tree(slong num_points);
	void clear_subproduct_tree(fq_poly_t **tree, slong num_points);
	
public:
	
	/**
	 * Evaluates {@code f} at the {@code n} elements {@code points} of the field
	 * {@code ctx}, and store the results in {@code results}. It is assumed that
	 * {@code n} < deg {@code f}. The algorithm used, is a divide and conquer 
	 * algorithm that builds a subproduct tree.   
	 * 
	 * @param results	the result of evaluation
	 * @param f			the given polynomial
	 * @param points	the given elements of the field {@code ctx}
	 * @param ctx		the given field
	 */
	void multipoint_eval(fq_t *results, const fq_poly_t f, const fq_t *points, slong num_points, const fq_ctx_t ctx);
	
};

#endif /* FQ_POLY_MULTIPPOINT_EVAL_H_ */
