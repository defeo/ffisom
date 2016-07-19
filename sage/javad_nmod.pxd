# distutils: libraries = javad_nmod
# distutils: library_dirs = .

from sage.libs.flint.nmod_poly cimport nmod_poly_t
from finite_field_flint_fq_nmod cimport FiniteField_flint_fq_nmod

cdef extern from "javad/nmod_poly_isom/ff_embedding.h":
    cdef cppclass FFEmbedding:
        FFEmbedding(nmod_poly_t, nmod_poly_t) except +
        void compute_generators(nmod_poly_t g1, nmod_poly_t g2, long lathr, long mpthr)
        void build_embedding(nmod_poly_t g1, nmod_poly_t g2)
        void compute_image(nmod_poly_t image, const nmod_poly_t f)
        void get_x_image(nmod_poly_t x)

cdef class FFEmbWrapper:
    cdef FFEmbedding *wrp
    cdef nmod_poly_t g1, g2, xim
    cdef FiniteField_flint_fq_nmod domain, codomain
    cdef int initialized
