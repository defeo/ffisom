# distutils: libraries = ellmul
# distutils: library_dirs = .

from sage.libs.flint.fmpz cimport fmpz_t
from sage.libs.flint.fq_nmod cimport fq_nmod_t, fq_nmod_ctx_t

cdef extern from "ellmul/fq_nmod_weierstrass_xz.h":
    ctypedef struct fq_nmod_weierstrass_xz_struct:
         pass
    ctypedef fq_nmod_weierstrass_xz_struct fq_nmod_weierstrass_xz_t[1]

    void fq_nmod_weierstrass_xz_set_ui(fq_nmod_weierstrass_xz_t, unsigned long, unsigned long, fq_nmod_ctx_t)
    void fq_nmod_weierstrass_xz_set_fq_nmod(fq_nmod_weierstrass_xz_t, fq_nmod_t, fq_nmod_t, fq_nmod_ctx_t)

    void fq_nmod_weierstrass_xz_init(fq_nmod_weierstrass_xz_t, fq_nmod_ctx_t)
    void fq_nmod_weierstrass_xz_clear(fq_nmod_weierstrass_xz_t, fq_nmod_ctx_t)

    void fq_nmod_weierstrass_xz_dbl(fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_weierstrass_xz_t, fq_nmod_ctx_t)
    void fq_nmod_weierstrass_xz_dadd(fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_weierstrass_xz_t, fq_nmod_ctx_t)
    void fq_nmod_weierstrass_xz_ladd(fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_weierstrass_xz_t, fq_nmod_ctx_t)

    void fq_nmod_weierstrass_xz_mul_ltr(fq_nmod_t, fq_nmod_t, fq_nmod_t, fq_nmod_t, fmpz_t, fq_nmod_weierstrass_xz_t, fq_nmod_ctx_t)

