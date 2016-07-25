/*
 * util.h
 *
 *  Created on: Feb 13, 2014
 *      Author: javad
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <flint/ulong_extras.h>

class Util {
    static const double EXT_COMP_EXPONENT = 0.2;

public:

    /**
     * Computes the order of {@code a} in the group $\mathbb{Z}/modulus\mathbb{Z}$.
     * It is assumed that {@code a} and {@code modulus} are coprime.
     */

    ulong compute_multiplicative_order(ulong a, ulong modulus);

    /**
     * Checks if an extension of degree {@degree} over $\mathbb{F}_p$ is small.
     * We consider a cyclotomic extension small if the multiplicative order $d$ of {@code p} in
     * the group $\mathbb{Z}/r\mathbb{Z}$ is such that $d < r^{2 - w}$ where $w$ is the
     * exponent of modular composition.
     */
    bool is_small_cyclotomic_ext(ulong degree, const ulong p);

};


#endif /* UTIL_H_ */
