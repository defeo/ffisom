/*
 * util.cpp
 *
 *  Created on: Feb 13, 2014
 *      Author: javad
 */



#include "util.h"
#include <math.h>

slong Util::compute_multiplicative_order(ulong a, ulong modulus) {

	n_factor_t factors;
	n_factor_init(&factors);

	ulong order = modulus - 1;
	n_factor(&factors, order, 1);

	slong temp = 1;
	for (slong i = 0; i < factors.num; i++) {
		while (temp == 1 && factors.exp[i] > 0) {
			order /= factors.p[i];
			factors.exp[i]--;
			temp = n_powmod(a, order, modulus);
		}

		if (temp != 1) {
			order *= factors.p[i];
			temp = 1;
		}
	}

	return order;
}


bool Util::is_small_cyclotomic_ext(ulong degree, ulong p){
	slong s = compute_multiplicative_order(p, degree);
	
	return (s < pow(degree, EXT_COMP_EXPONENT));
}
