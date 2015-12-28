#clib pari gmp
#from sage.env import SAGE_ENV
#include SAGE_ENV["SITE_PACKAGES"]+"sage/ext/interrupt.pxi"
include "sage/ext/interrupt.pxi"

from sage.libs.gmp.mpz cimport mpz_fits_ulong_p, mpz_get_ui

from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer

import sage.libs.pari.pari_instance
from sage.libs.pari.pari_instance cimport PariInstance
from sage.libs.pari.types cimport GEN, t_VECSMALL, evalvarn, lg
from sage.libs.pari.paridecl cimport avma, cgetg, Flx_ffisom, Flx_ffintersect, ZX_to_Flx
from sage.libs.pari.paridecl cimport pari_printf, Flx_is_irred

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
    cdef PariInstance PI = sage.libs.pari.pari_instance.pari
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
    sig_off()
    sig_on()
    a = [(<unsigned long *>SP)[i] for i in xrange(2, lg(SP))]
    b = [(<unsigned long *>SQ)[i] for i in xrange(2, lg(SQ))]
    PI.clear_stack()

    return domain(a), codomain(b)

def find_emb(domain, codomain):
    cdef PariInstance PI = sage.libs.pari.pari_instance.pari
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

    P = intlist_to_Flx(domain.modulus().change_ring(ZZ).list())
    Q = intlist_to_Flx(codomain.modulus().change_ring(ZZ).list())
    pari_printf("%Ps\n",P);
    pari_printf("%Ps\n",Q);
    sig_on()
    xim = Flx_ffisom(P, Q, l)
    a = [(<unsigned long *>xim)[i] for i in xrange(2, lg(xim))]
    sig_off()
    sig_on()
    PI.clear_stack()

    return codomain(a)
