#clib pari gmp
include "cysignals/signals.pxi" 

from sage.libs.gmp.mpz cimport mpz_fits_ulong_p, mpz_get_ui

from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer

from sage.libs.pari.stack cimport clear_stack
from sage.libs.pari.types cimport GEN, t_VECSMALL, evalvarn, lg
from sage.libs.pari.paridecl cimport avma, cgetg, gel, Flx_ffisom, Flx_ffintersect, ZX_to_Flx, Flx_factorff_irred

# Dirty stack
cdef GEN intlist_to_Flx(list L):
    cdef int i, l
    cdef GEN a
    l = len(L) + 2
    a = cgetg(l, t_VECSMALL);
    a[1] = evalvarn(0)
    for i in xrange(2, l):
        a[i] = L[i-2]
        #a[i] = mpz_get_ui((<Integer> L[i-2]).value)
    return a

def find_gen(domain, r = 0):
    raise NotImplementedError

def find_gens(domain, codomain, int r = 0):
    cdef int i
    #cdef Integer p
    cdef unsigned long l
    cdef GEN P, Q
    cdef GEN SP, SQ
    cdef list a, b


    #p = domain.characteristic()
    #if not mpz_fits_ulong_p(p.value):
    #    raise ValueError
    #l = mpz_get_ui(p.value)
    l = domain.characteristic()

    if r == 0:
        r = domain.degree()
    #    r = mpz_get_ui((<Integer>domain.degree()).value)

    sig_on()
    P = intlist_to_Flx(domain.modulus().change_ring(ZZ).list())
    Q = intlist_to_Flx(codomain.modulus().change_ring(ZZ).list())
    Flx_ffintersect(P, Q, r, l, &SP, &SQ, NULL, NULL)
    a = [(<unsigned long *> SP)[i] for i in xrange(2, lg(SP))]
    b = [(<unsigned long *> SQ)[i] for i in xrange(2, lg(SQ))]
    clear_stack()

    return domain(a), codomain(b)

def find_emb(domain, codomain):
    cdef int i
    #cdef Integer p
    cdef unsigned long l
    cdef GEN P, Q
    cdef GEN xim
    cdef list a

    #p = domain.characteristic()
    #if not mpz_fits_ulong_p(p.value):
    #    raise ValueError
    #l = mpz_get_ui(p.value)
    l = domain.characteristic()

    sig_on()
    P = intlist_to_Flx(domain.modulus().change_ring(ZZ).list())
    Q = intlist_to_Flx(codomain.modulus().change_ring(ZZ).list())
    xim = Flx_ffisom(P, Q, l)
    a = [(<unsigned long *> xim)[i] for i in xrange(2, lg(xim))]
    clear_stack()

    return codomain(a)

def find_root(poly, domain):
    cdef int i
    #cdef Integer p
    cdef unsigned long l
    cdef GEN P, Q
    cdef GEN factors, root
    cdef list a

    #p = domain.characteristic()
    #if not mpz_fits_ulong_p(p.value):
    #    raise ValueError
    #l = mpz_get_ui(p.value)
    l = domain.characteristic()

    sig_on()
    P = intlist_to_Flx(poly.change_ring(ZZ).list())
    Q = intlist_to_Flx(domain.modulus().change_ring(ZZ).list())
    factors = Flx_factorff_irred(P, Q, l)
    root = <GEN> (<GEN> factors[2])[2]
    a = [l-(<unsigned long *> root)[i] for i in xrange(2, lg(root))]
    clear_stack()

    return domain(a)
