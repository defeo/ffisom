from sage.libs.flint.types cimport fq_nmod_ctx_struct

from sage.rings.finite_rings.finite_field_base cimport FiniteField

cdef class FiniteField_flint_fq_nmod(FiniteField):
    cdef fq_nmod_ctx_struct *_ctx
    #cdef fq_nmod_ctx_t _ctx
    cdef int _ctx_initialized

    cdef public object _modulus
    cdef public object _degree
    cdef public object _gen
    cdef public int __hash
