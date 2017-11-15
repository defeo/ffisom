from cysignals.signals cimport sig_on, sig_off

from sage.libs.flint.types cimport fq_nmod_t, fq_nmod_ctx_struct
from sage.libs.flint.fq_nmod cimport fq_nmod_init, fq_nmod_clear, fq_nmod_set_si
from finite_field_flint_fq_nmod cimport FiniteField_flint_fq_nmod
from element_flint_fq_nmod cimport FiniteFieldElement_flint_fq_nmod
 
def find_root(poly, FiniteField_flint_fq_nmod k):
    cdef fq_nmod_ctx_struct *ctx
    cdef fq_nmod_t coeff, root
    cdef fq_nmod_poly_t f, factor
    cdef flint_rand_t state
    cdef int found
    cdef FiniteFieldElement_flint_fq_nmod a

    ctx = k._ctx

    fq_nmod_init(coeff, ctx)
    fq_nmod_poly_init(f, ctx)
    for i, c in enumerate(poly.list()):
        fq_nmod_set_si(coeff, <long> c, ctx)
        fq_nmod_poly_set_coeff(f, <long> i, coeff, ctx)
    fq_nmod_clear(coeff, ctx)

    fq_nmod_poly_init(factor, ctx)
    found = 0
    flint_randinit(state)
    sig_on()
    while fq_nmod_poly_degree(factor, ctx) != 1:
        while not found:
            found = fq_nmod_poly_factor_equal_deg_prob(factor, state, f, 1, ctx)
        fq_nmod_poly_set(f, factor, ctx)
        found = 0
    sig_off()
    flint_randclear(state)

    fq_nmod_init(root, ctx)
    fq_nmod_poly_get_coeff(root, factor, 0, ctx)
    a = k(0)
    a.set_from_fq_nmod(root)
    a = -a
    fq_nmod_clear(root, ctx)

    fq_nmod_poly_clear(factor, ctx)
    fq_nmod_poly_clear(f, ctx)

    return a
