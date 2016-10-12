# distutils: include_dirs = kummer_c++_flint
# distutils: libraries = kummer
# distutils: library_dirs = .

from sage.libs.flint.nmod_poly cimport nmod_poly_t
from finite_field_flint_fq_nmod cimport FiniteField_flint_fq_nmod

cdef extern from "kummer_c++_flint/ff_embedding.h":
    cdef cppclass FFEmbedding:
        FFEmbedding(nmod_poly_t, nmod_poly_t, long) except +
        void compute_generators(nmod_poly_t g1, nmod_poly_t g2)
        void build_embedding(nmod_poly_t g1, nmod_poly_t g2)
        void compute_image(nmod_poly_t image, const nmod_poly_t f)
        void get_x_image(nmod_poly_t x)

cdef extern from "kummer_c++_flint/ff_isom_prime_power_ext.h":
    cdef enum:
        FORCE_LINALG
        FORCE_MODCOMP
        FORCE_COFACTOR
        FORCE_MPE
        FORCE_NONE

cdef class FFEmbWrapper:
    cdef FFEmbedding *wrp
    cdef nmod_poly_t g1, g2, xim
    cdef FiniteField_flint_fq_nmod domain, codomain
    cdef int initialized
