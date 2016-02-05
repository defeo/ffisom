#distutils: libraries = flint

from sage.libs.flint.types cimport slong#, flint_rand_t
from sage.libs.flint.types cimport fq_nmod_ctx_t, fq_nmod_t

cdef extern from "flint/flint.h":
    ctypedef struct flint_rand_s:
        pass

    ctypedef flint_rand_s flint_rand_t[1]

    void flint_randinit(flint_rand_t state)
    void flint_randclear(flint_rand_t state)

cdef extern from "flint/fq_nmod_poly.h":
    ctypedef struct fq_nmod_poly_struct:
        pass

    ctypedef fq_nmod_poly_struct fq_nmod_poly_t[1]

    void fq_nmod_poly_init(fq_nmod_poly_t poly, fq_nmod_ctx_t ctx)
    void fq_nmod_poly_fit_length(fq_nmod_poly_t poly, slong len,
                              const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_set_length(fq_nmod_poly_t poly, slong len,
                              const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_clear(fq_nmod_poly_t poly, fq_nmod_ctx_t ctx)
    void fq_nmod_poly_normalise(fq_nmod_poly_t, fq_nmod_ctx_t)

    long fq_nmod_poly_degree(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)

    void fq_nmod_poly_set(fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2,
                      const fq_nmod_ctx_t ctx)

    void fq_nmod_poly_get_coeff(fq_nmod_t x, const fq_nmod_poly_t poly, slong n,
                            const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_set_coeff(fq_nmod_poly_t poly, slong n, const fq_nmod_t x,
                            const fq_nmod_ctx_t ctx)

    void fq_nmod_poly_zero(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
    void fq_nmod_poly_one(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)

    int fq_nmod_poly_print_pretty(const fq_nmod_poly_t poly, const char *x,
                         const fq_nmod_ctx_t ctx)

cdef extern from "flint/fq_nmod_poly_factor.h":
    int fq_nmod_poly_factor_equal_deg_prob(fq_nmod_poly_t factor, flint_rand_t state,
                                  const fq_nmod_poly_t pol, slong d,
                                  const fq_nmod_ctx_t ctx)
