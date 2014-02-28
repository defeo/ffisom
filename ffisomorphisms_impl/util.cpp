/*
 * util.cpp
 *
 *  Created on: Feb 13, 2014
 *      Author: javad
 */



#include "util.h"
#include <math.h>

slong Util::compute_multiplicative_order(slong a, slong modulus) {

	n_factor_t factors;
	n_factor_init(&factors);

	slong order = modulus - 1;
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

slong Util::compute_multiplicative_order(const fmpz_t a, slong modulus){
	return compute_multiplicative_order(fmpz_fdiv_ui(a, modulus), modulus);
}

bool Util::is_small_cyclotomic_ext(slong degree, const fmpz_t p){
	slong s = compute_multiplicative_order(p, degree);
	
	return (s < pow(degree, EXT_COMP_EXPONENT));
}
