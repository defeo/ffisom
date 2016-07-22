from sage.libs.flint.types cimport fq_nmod_t, fq_nmod_ctx_struct

from sage.rings.finite_rings.element_base cimport FinitePolyExtElement

cdef class FiniteFieldElement_flint_fq_nmod(FinitePolyExtElement):
    cdef fq_nmod_t val
    cdef int initialized
    cdef fq_nmod_ctx_struct *_cparent
    cdef FiniteFieldElement_flint_fq_nmod _new(FiniteFieldElement_flint_fq_nmod self)
    cdef void set_from_fq_nmod(FiniteFieldElement_flint_fq_nmod self, fq_nmod_t val) except *
    cdef void construct_from(FiniteFieldElement_flint_fq_nmod self, object x) except *
