// copyright blablabla GPL v3 or later at your option

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mod_poly.h"
#include "flint/fq.h"

// support nails?
// make the curve types ctx like? and add a point type?
// store finit field ctx within curve ctx?
// make the multiplication routine generic? using function pointers?
// pass tmp space to functions?
// use windows?
// return z==0 for bdl, dadd, ladd, mul functions?
// add z=1 formulas?
// add conversion function, extract y from x?

// Montgomery XZ type
typedef struct {
    fq_t a;
    fq_t a24;
    fq_t b;
//    fq_ctx_t K;
} montgomery_xz_struct;

typedef montgomery_xz_struct montgomery_xz_t[1];

void montgomery_xz_init(montgomery_xz_t E, fq_ctx_t K)
{
    fq_init(E->a, K);
    fq_init(E->a24, K);
    fq_init(E->b, K);
}

// requires the inverse of 4 in K
void montgomery_xz_set_ui(montgomery_xz_t E, ulong a, ulong b, fq_ctx_t K, fq_t inv4)
{
    fq_set_ui(E->a, a, K);
    fq_set_ui(E->a24, 2, K);
    fq_add(E->a24, E->a, E->a24, K);
    fq_mul(E->a24, E->a24, inv4, K);
    fq_set_ui(E->b, b, K);
}

void montgomery_xz_clear(montgomery_xz_t E, fq_ctx_t K)
{
    fq_clear(E->a, K);
    fq_clear(E->a24, K);
    fq_clear(E->b, K);
}

// Montgomery XZ doubling
// P1 |-> P3=[2]P1
// dbl-1987-m-3
void montgomery_xz_dbl(fq_t x3, fq_t z3,
                       const fq_t x1, const fq_t z1,
                       const montgomery_xz_t E,
                       const fq_ctx_t K)
{
    fq_t a, b, aa, bb;
    fq_init(a, K);
    fq_init(b, K);
    fq_init(aa, K);
    fq_init(bb, K);

    // A = X1+Z1
    fq_add(a, x1, z1, K);
    // AA = A²
    fq_sqr(aa, a, K);
    // B = X1-Z1
    fq_sub(b, x1, z1, K);
    // BB = B²
    fq_sqr(bb, b, K);
    // C = AA-BB
    fq_sub(x3, aa, bb, K);
    // Z3 = C*(BB+a24*C)
    fq_mul(a, E->a24, x3, K);
    fq_add(b, bb, a, K);
    fq_mul(z3, x3, b, K);
    // X3 = AA*BB
    fq_mul(x3, aa, bb, K);

    fq_clear(a, K);
    fq_clear(b, K);
    fq_clear(aa, K);
    fq_clear(bb, K);
}

// Montgomery XZ differential addition
// P1=P3-P2, P2, P3 |-> P5=P2+P3
// dadd-1987-m-3
void montgomery_xz_dadd(fq_t x5, fq_t z5,
                        const fq_t x3, const fq_t z3,
                        const fq_t x2, const fq_t z2,
                        const fq_t x1, const fq_t z1,
                        const montgomery_xz_t E,
                        const fq_ctx_t K)
{
    fq_t a, b, c, d, da, cb;
    fq_init(a, K);
    fq_init(b, K);
    fq_init(c, K);
    fq_init(d, K);
    fq_init(da, K);
    fq_init(cb, K);

    // A = X2+Z2
    fq_add(a, x2, z2, K);
    // B = X2-Z2
    fq_sub(b, x2, z2, K);
    // C = X3+Z3
    fq_add(c, x3, z3, K);
    // D = X3-Z3
    fq_sub(d, x3, z3, K);
    // DA = D*A
    fq_mul(da, d, a, K);
    // CB = C*B
    fq_mul(cb, c, b, K);
    // X5 = Z1*(DA+CB)²
    fq_add(x5, da, cb, K);
    fq_sqr(a, x5, K);
    fq_mul(x5, z1, a, K);
    // Z5 = X1*(DA-CB)²
    fq_sub(z5, da, cb, K);
    fq_sqr(a, z5, K);
    fq_mul(z5, x1, a, K);

    fq_clear(a, K);
    fq_clear(b, K);
    fq_clear(c, K);
    fq_clear(d, K);
    fq_clear(da, K);
    fq_clear(cb, K);
}

// Montgomery XZ ladder
// P1=P3-P2, P2, P3 |-> P4=[2]P2, P5=P2+P3
// ladd-1987-m-3
void montgomery_xz_ladd(fq_t x5, fq_t z5, fq_t x4, fq_t z4,
                        const fq_t x3, const fq_t z3,
                        const fq_t x2, const fq_t z2,
                        const fq_t x1, const fq_t z1,
                        const montgomery_xz_t E,
                        const fq_ctx_t K)
{
    fq_t a, b, c, d, e, aa, bb;
    fq_init(a, K);
    fq_init(b, K);
    fq_init(c, K);
    fq_init(d, K);
    fq_init(e, K);
    fq_init(aa, K);
    fq_init(bb, K);

    // A = X2+Z2
    fq_add(a, x2, z2, K);
    // AA = A²
    fq_sqr(aa, a, K);
    // B = X2-Z2
    fq_sub(b, x2, z2, K);
    // BB = B²
    fq_sqr(bb, b, K);
    // E = AA-BB
    fq_sub(e, aa, bb, K);
    // X4 = AA*BB
    fq_mul(x4, aa, bb, K);
    // Z4 = E*(BB+a24*E)
    fq_mul(z4, E->a24, e, K);
    fq_add(aa, bb, z4, K);
    fq_mul(z4, e, aa, K);
    // C = X3+Z3
    fq_add(c, x3, z3, K);
    // D = X3-Z3
    fq_sub(d, x3, z3, K);
    // DA = D*A
    fq_mul(aa, d, a, K);
    // CB =  C*B
    fq_mul(bb, c, b, K);
    // X5 = Z1*(DA+CB)²
    fq_add(x5, aa, bb, K);
    fq_sqr(e, x5, K);
    fq_mul(x5, z1, e, K);
    // Z5 = X1*(DA-CB)²
    fq_sub(z5, aa, bb, K);
    fq_sqr(e, z5, K);
    fq_mul(z5, x1, e, K);

    fq_clear(a, K);
    fq_clear(b, K);
    fq_clear(c, K);
    fq_clear(d, K);
    fq_clear(e, K);
    fq_clear(aa, K);
    fq_clear(bb, K);
}

// Montgomery XZ point multiplication
// Left-to-right binary exponentiation
montgomery_xz_mul_ltr(fq_t x, fq_t z,
                       const fq_t x1, const fq_t z1,
                       const fmpz_t m,
                       const montgomery_xz_t E,
                       const fq_ctx_t K)
{
    if (fmpz_is_zero(m))
    {
        fq_one(x, K);
        fq_zero(z, K);
    }
    else if (fmpz_is_one(m))
    {
        fq_set(x, x1, K);
        fq_set(z, z1, K);
    }
    else
    {
        ulong bit;
        fq_t x2, z2, x3, z3, x4, z4, x5, z5;


        fq_init(x2, K);
        fq_init(z2, K);
        fq_init(x3, K);
        fq_init(z3, K);
        fq_init(x4, K);
        fq_init(z4, K);
        fq_init(x5, K);
        fq_init(z5, K);

        // P
        fq_set(x2, x1, K);
        fq_set(z2, z1, K);
        // 2*P
        montgomery_xz_dbl(x3, z3, x1, z1, E, K);

        for (bit = fmpz_bits(m) - 2; bit > 0; bit--)
        {
            if (fmpz_tstbit(m, bit))
            {
                // P, [k+1]P, [k]P -> [2k+2]P, [2k+1]P
                montgomery_xz_ladd(x5, z5, x4, z4, x2, z2, x3, z3, x1, z1, E, K);
                fq_swap(x2, x5, K);
                fq_swap(z2, z5, K);
                fq_swap(x3, x4, K);
                fq_swap(z3, z4, K);
            }
            else
            {
                // P, [k]P, [k+1]P -> [2k]P, [2k+1]P
                montgomery_xz_ladd(x5, z5, x4, z4, x3, z3, x2, z2, x1, z1, E, K);
                fq_swap(x2, x4, K);
                fq_swap(z2, z4, K);
                fq_swap(x3, x5, K);
                fq_swap(z3, z5, K);
            }
        }

        // Last iteration, no need for a ladder
        {
            if (fmpz_tstbit(m, 0))
            {
                montgomery_xz_dadd(x, z, x3, z3, x2, z2, x1, z1, E, K);
            }
            else
            {
                montgomery_xz_dbl(x, z, x2, z2, E, K);
            }
        }

        fq_clear(x2, K);
        fq_clear(z2, K);
        fq_clear(x3, K);
        fq_clear(z3, K);
        fq_clear(x4, K);
        fq_clear(z4, K);
        fq_clear(x5, K);
        fq_clear(z5, K);
    }
}

// Weierstrass XZ type
typedef struct {
    fq_t a;
    fq_t b;
    fq_t b2;
    fq_t b4;
//    fq_ctx_t K;
} weierstrass_xz_struct;

typedef weierstrass_xz_struct weierstrass_xz_t[1];

void weierstrass_xz_init(weierstrass_xz_t E, fq_ctx_t K)
{
    fq_init(E->a, K);
    fq_init(E->b, K);
    fq_init(E->b2, K);
    fq_init(E->b4, K);
}

void weierstrass_xz_set_ui(weierstrass_xz_t E, ulong a, ulong b, fq_ctx_t K)
{
    fq_set_ui(E->a, a, K);
    fq_set_ui(E->b, b, K);
    fq_mul_ui(E->b2, E->b, 2, K);
    fq_mul_ui(E->b4, E->b, 4, K);
}

void weierstrass_xz_clear(weierstrass_xz_t E, fq_ctx_t K)
{
    fq_clear(E->a, K);
    fq_clear(E->b, K);
    fq_clear(E->b2, K);
    fq_clear(E->b4, K);
}

// Weierstrass XZ doubling
// P1 |-> P3=[2]P1
// dbl-2002-bj-3
void weierstrass_xz_dbl(fq_t x3, fq_t z3,
                       const fq_t x1, const fq_t z1,
                       const weierstrass_xz_t E,
                       const fq_ctx_t K)
{
    fq_t xx, zz, a, azz;
    fq_init(xx, K);
    fq_init(zz, K);
    fq_init(a, K);
    fq_init(azz, K);

    // XX = X1²
    fq_sqr(xx, x1, K);
    // ZZ = Z1²
    fq_sqr(zz, z1, K);
    // A = 2*((X1+Z1)²-XX-ZZ)
    fq_add(azz, x1, z1, K);
    fq_sqr(a, azz, K);
    fq_sub(a, a, xx, K);
    fq_sub(a, a, zz, K);
    fq_mul_ui(a, a, 2, K);
    // aZZ = a*ZZ
    fq_mul(azz, E->a, zz, K);
    // X3 = (XX-aZZ)²-b2*A*ZZ
    fq_sub(x3, xx, azz, K);
    fq_sqr(x3, x3, K);
    fq_mul(z3, a, zz, K);
    fq_mul(z3, E->b2, z3, K);
    fq_sub(x3, x3, z3, K);
    // Z3 = A*(XX+aZZ)+b4*ZZ²
    fq_add(xx, xx, azz, K);
    fq_mul(azz, a, xx, K);
    fq_sqr(a, zz, K);
    fq_mul(zz, E->b4, a, K);
    fq_add(z3, azz, zz, K);

    fq_clear(xx, K);
    fq_clear(zz, K);
    fq_clear(a, K);
    fq_clear(azz, K);
}

// Weierstrass XZ differential addition
// P1=P3-P2, P2, P3 |-> P5=P2+P3
// dadd-2002-it-3
void weierstrass_xz_dadd(fq_t x5, fq_t z5,
                         const fq_t x3, const fq_t z3,
                         const fq_t x2, const fq_t z2,
                         const fq_t x1, const fq_t z1,
                         const weierstrass_xz_t E,
                         const fq_ctx_t K)
{
    fq_t a, b, c, d, e, f;
    fq_init(a, K);
    fq_init(b, K);
    fq_init(c, K);
    fq_init(d, K);
    fq_init(e, K);
    fq_init(f, K);

    // A = X2*X3
    fq_mul(a, x2, x3, K);
    // B = Z2*Z3
    fq_mul(b, z2, z3, K);
    // C = X2*Z3
    fq_mul(c, x2, z3, K);
    // D = X3*Z2
    fq_mul(d, x3, z2, K);
    // E = (A-a*B)²
    fq_mul(e, E->a, b, K);
    fq_sub(f, a, e, K);
    fq_sqr(e, f, K);
    // F = 4*b*B*(C+D)
    fq_add(f, c, d, K);
    fq_mul(a, b, f, K);
    fq_mul(b, E->b, a, K);
    fq_mul_ui(f, b, 4, K);
    // X5 = Z1*(E-F)
    fq_sub(a, e, f, K);
    fq_mul(x5, z1, a, K);
    // Z5 = X1*(C-D)²
    fq_sub(a, c, d, K);
    fq_sqr(b, a, K);
    fq_mul(z5, x1, b, K);

    fq_clear(a, K);
    fq_clear(b, K);
    fq_clear(c, K);
    fq_clear(d, K);
    fq_clear(e, K);
    fq_clear(f, K);
}

// Weierstrass XZ ladder
// P1=P3-P2, P2, P3 |-> P4=[2]P2, P5=P2+P3
// ladd-2002-it-3
void weierstrass_xz_ladd(fq_t x5, fq_t z5, fq_t x4, fq_t z4,
                         const fq_t x3, const fq_t z3,
                         const fq_t x2, const fq_t z2,
                         const fq_t x1, const fq_t z1,
                         const weierstrass_xz_t E,
                         const fq_ctx_t K)
{
    fq_t xx, zz, e, azz, a, b, c, d;
    fq_init(xx, K);
    fq_init(zz, K);
    fq_init(azz, K);
    fq_init(e, K);
    fq_init(a, K);
    fq_init(b, K);
    fq_init(c, K);
    fq_init(d, K);

    // XX = X2²
    fq_sqr(xx, x2, K);
    // ZZ = Z2²
    fq_sqr(zz, z2, K);
    // aZZ = a*ZZ
    fq_mul(azz, E->a, zz, K);
    // E = (X2+Z2)²-XX-ZZ
    fq_add(a, x2, z2, K);
    fq_sqr(b, a, K);
    fq_add(a, xx, zz, K);
    fq_sub(e, b, a, K);
    // X4 = (XX-aZZ)²-b4*E*ZZ
    fq_sub(a, xx, azz, K);
    fq_sqr(b, a, K);
    fq_mul(a, e, zz, K);
    fq_mul(c, a, E->b4, K);
    fq_sub(x4, b, c, K);
    // Z4 = 2*E*(XX+aZZ)+b4*ZZ²
    fq_add(a, xx, azz, K);
    fq_mul(b, a, e, K);
    fq_mul_ui(c, b, 2, K);
    fq_sqr(a, zz, K);
    fq_mul(b, E->b4, a, K);
    fq_add(z4, c, b, K);
    // A = X2*X3
    fq_mul(a, x2, x3, K);
    // B = Z2*Z3
    fq_mul(b, z2, z3, K);
    // C = X2*Z3
    fq_mul(c, x2, z3, K);
    // D = X3*Z2
    fq_mul(d, x3, z2, K);
    // X5 = Z1*((A-a*B)²-b4*B*(C+D))
    fq_mul(azz, E->a, b, K);
    fq_sub(x5, a, azz, K);
    fq_sqr(a, x5, K);
    fq_add(x5, c, d, K);
    fq_mul(e, b, x5, K);
    fq_mul(b, E->b4, e, K);
    fq_sub(e, a, b, K);
    fq_mul(x5, z1, e, K);
    // Z5 = X1*(C-D)²
    fq_sub(a, c, d, K);
    fq_sqr(b, a, K);
    fq_mul(z5, x1, b, K);

    fq_clear(xx, K);
    fq_clear(zz, K);
    fq_clear(azz, K);
    fq_clear(e, K);
    fq_clear(a, K);
    fq_clear(b, K);
    fq_clear(c, K);
    fq_clear(d, K);
}

// Make the following generic? using hackish pointers?
// Weierstrass XZ point multiplication
// Left-to-right binary exponentiation
weierstrass_xz_mul_ltr(fq_t x, fq_t z,
                       const fq_t x1, const fq_t z1,
                       const fmpz_t m,
                       const weierstrass_xz_t E,
                       const fq_ctx_t K)
{
    if (fmpz_is_zero(m))
    {
        fq_one(x, K);
        fq_zero(z, K);
    }
    else if (fmpz_is_one(m))
    {
        fq_set(x, x1, K);
        fq_set(z, z1, K);
    }
    else
    {
        ulong bit;
        fq_t x2, z2, x3, z3, x4, z4, x5, z5;


        fq_init(x2, K);
        fq_init(z2, K);
        fq_init(x3, K);
        fq_init(z3, K);
        fq_init(x4, K);
        fq_init(z4, K);
        fq_init(x5, K);
        fq_init(z5, K);

        // P
        fq_set(x2, x1, K);
        fq_set(z2, z1, K);
        // 2*P
        weierstrass_xz_dbl(x3, z3, x1, z1, E, K);

        for (bit = fmpz_bits(m) - 2; bit > 0; bit--)
        {
            if (fmpz_tstbit(m, bit))
            {
                // P, [k+1]P, [k]P -> [2k+2]P, [2k+1]P
                weierstrass_xz_ladd(x5, z5, x4, z4, x2, z2, x3, z3, x1, z1, E, K);
                fq_swap(x2, x5, K);
                fq_swap(z2, z5, K);
                fq_swap(x3, x4, K);
                fq_swap(z3, z4, K);
            }
            else
            {
                // P, [k]P, [k+1]P -> [2k]P, [2k+1]P
                weierstrass_xz_ladd(x5, z5, x4, z4, x3, z3, x2, z2, x1, z1, E, K);
                fq_swap(x2, x4, K);
                fq_swap(z2, z4, K);
                fq_swap(x3, x5, K);
                fq_swap(z3, z5, K);
            }
        }

        // Last iteration, no need for a ladder
        {
            if (fmpz_tstbit(m, 0))
            {
                weierstrass_xz_dadd(x, z, x3, z3, x2, z2, x1, z1, E, K);
            }
            else
            {
                weierstrass_xz_dbl(x, z, x2, z2, E, K);
            }
        }

        fq_clear(x2, K);
        fq_clear(z2, K);
        fq_clear(x3, K);
        fq_clear(z3, K);
        fq_clear(x4, K);
        fq_clear(z4, K);
        fq_clear(x5, K);
        fq_clear(z5, K);
    }
}

// Test routine
int main(int argc, char* argv[])
{
    unsigned long n;
    n = 300;

    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, 1031);

    fmpz_t q;
    fmpz_init(q);
    fmpz_pow_ui(q, p, n);
    //fmpz_set_ui(q, 4);
    fmpz_print(q);
    printf("\n");

    fmpz_mod_poly_t f;
    fmpz_mod_poly_init(f, p);
    fmpz_mod_poly_set_coeff_ui(f, 300, 1);
    fmpz_mod_poly_set_coeff_ui(f, 1, 1);
    fmpz_mod_poly_set_coeff_ui(f, 0, 321);

    fq_ctx_t K;
    fq_ctx_init_modulus(K, f, "a");

    fq_t inv4;
    fq_init(inv4, K);
    fq_set_ui(inv4, 4, K);
    fq_inv(inv4, inv4, K);

    fq_t x1, z1;
    fq_init(x1, K);
    fq_init(z1, K);
    fq_one(z1, K);
    fq_gen(x1, K);
    fq_add(x1, x1, z1, K);

    fq_t x, z;
    fq_init(x, K);
    fq_init(z, K);

    clock_t t;

    // Weierstrass XZ
    weierstrass_xz_t W;
    weierstrass_xz_init(W, K);
    weierstrass_xz_set_ui(W, 207, 693, K);

    t = clock();
    weierstrass_xz_mul_ltr(x, z, x1, z1, q, W, K);
    printf("Weierstrass XZ: %lf\n", ((double) clock() - t) / CLOCKS_PER_SEC);
    fq_inv(z1, z, K);
    fq_mul(x1, x, z1, K);
    fq_print_pretty(x1, K);
    printf("\n");

    weierstrass_xz_clear(W, K);

    // Montgomery XZ
    montgomery_xz_t M;
    montgomery_xz_init(M, K);
    montgomery_xz_set_ui(M, 207, 693, K, inv4);

    t = clock();
    montgomery_xz_mul_ltr(x, z, x1, z1, q, M, K);
    printf("Montgomery XZ: %lf\n", ((double) clock() - t) / CLOCKS_PER_SEC);
    fq_inv(z1, z, K);
    fq_mul(x1, x, z1, K);
    fq_print_pretty(x1, K);
    printf("\n");

    montgomery_xz_clear(M, K);

    fq_clear(x, K);
    fq_clear(z, K);
    fq_clear(x1, K);
    fq_clear(z1, K);
    fq_ctx_clear(K);
    fmpz_mod_poly_clear(f);
    fmpz_clear(q);
    fmpz_clear(p);

    return 0;
}
