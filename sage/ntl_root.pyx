# distutils: language = c++

include "sage/ext/interrupt.pxi"

def find_root(poly, k):
    cdef zz_pEPush pEPushed
    cdef zz_pPush pPushed

    p = k.characteristic()
    zz_p_init(p)
    cdef zz_pX modulus
    for i, c in enumerate(k.modulus().list()):
        SetCoeff(modulus, <long> i, <long> c)

    zz_pE_init(modulus)
    cdef zz_pEX f
    for i, c in enumerate(poly.list()):
        SetCoeff(f, <long> i, <long> c)

    sig_on()
    cdef zz_pE root = FindRoot(f)
    sig_off()

    cdef zz_pX root_poly = rep(root)

    return k([rep(coeff(root_poly, i)) for i in xrange(deg(root_poly) + 1)])
