# distutils: language = c++
include "sage/ext/stdsage.pxi"
include "cysignals/signals.pxi"

from finite_field_flint_fq_nmod cimport FiniteField_flint_fq_nmod
from element_flint_fq_nmod cimport FiniteFieldElement_flint_fq_nmod
from sage.libs.flint.nmod_poly cimport *
from sage.libs.flint.fq_nmod cimport *

cdef class FFEmbWrapper:
    def __cinit__(self, FiniteField_flint_fq_nmod k1, FiniteField_flint_fq_nmod k2):
        self.wrp = new FFEmbedding(k1._ctx.modulus, k2._ctx.modulus)
        nmod_poly_init(self.g1, k1._ctx.modulus.mod.n)
        nmod_poly_init(self.g2, k2._ctx.modulus.mod.n)
        nmod_poly_init(self.xim, k2._ctx.modulus.mod.n)

    def __init__(self, FiniteField_flint_fq_nmod k1, FiniteField_flint_fq_nmod k2):
        self.domain = k1
        self.codomain = k2
        self.initialized = 0

    def __dealloc__(self):
        del self.wrp
        nmod_poly_clear(self.g1)
        nmod_poly_clear(self.g2)
        nmod_poly_clear(self.xim)

    def compute_gens(self):
        self.wrp.compute_generators(self.g1, self.g2)
        self.initialized = 1

    def get_gens(self):
        cdef FiniteFieldElement_flint_fq_nmod g1, g2

        if self.initialized < 1:
            self.compute_gens()

        g1 = self.domain(0)
        g1.set_from_fq_nmod(<fq_nmod_t>(self.g1))
        g2 = self.codomain(0)
        g2.set_from_fq_nmod(<fq_nmod_t>(self.g2))

        return (g1, g2)

    def compute_emb(self):
        if self.initialized < 1:
            self.compute_gens()

        self.wrp.build_embedding(self.g1, self.g2)
        self.wrp.get_x_image(self.xim)
        self.initialized = 2

    def get_emb(self):
        cdef FiniteFieldElement_flint_fq_nmod xim

        if self.initialized < 2:
            self.compute_emb()

        xim = self.codomain(0)
        xim.set_from_fq_nmod(<fq_nmod_t>(self.xim))

        return xim

def find_gen(k1, r = 0):
    raise NotImplementedError

def find_gens(k1, k2, int r = 0):
    return FFEmbWrapper(k1, k2).get_gens()

def find_emb(k1, k2):
    return FFEmbWrapper(k1, k2).get_emb()
