/*
 * fq_nmod_poly_multippoint_eval.h
 *
 *  Created on: Dec 5, 2013
 *      Author: javad
 */

#ifndef fq_nmod_POLY_MULTIPPOINT_EVAL_H_
#define fq_nmod_POLY_MULTIPPOINT_EVAL_H_

#include <flint/fq_nmod_poly.h>


/**
 * This class is fast multipoint evaluation over finite fields.
 */
class fq_nmodPolyEval {

	const fq_nmod_ctx_struct *ctx;

	void build_subproduct_tree(fq_nmod_poly_t **tree, const fq_nmod_t *points, slong num_points);
	void go_down_subproduct_tree(fq_nmod_t *results, const fq_nmod_poly_t f, fq_nmod_poly_t **tree, slong num_points);
	fq_nmod_poly_t ** init_subproduct_tree(slong num_points);
	void clear_subproduct_tree(fq_nmod_poly_t **tree, slong num_points);
	
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
	void multipoint_eval(fq_nmod_t *results, const fq_nmod_poly_t f, const fq_nmod_t *points, slong num_points, const fq_nmod_ctx_t ctx);
	
};

#endif /* fq_nmod_POLY_MULTIPPOINT_EVAL_H_ */
