# distutils: libraries = ntl

cdef extern from "NTL/lzz_p.h" namespace "NTL":
    cdef cppclass zz_p:
        zz_p()
        zz_p(long a)

    void zz_p_init "NTL::zz_p::init"(long p)

    long rep(zz_p a)

    cdef cppclass zz_pPush:
        pass

cdef extern from "NTL/lzz_pX.h" namespace "NTL":
    cdef cppclass zz_pX:
        pass

    long deg(const zz_pX& a)
    const zz_p coeff(const zz_pX& a, long i)
    void SetCoeff(zz_pX& x, long i, long a)
    void SetCoeff(zz_pX& x, long i)
    void SetX(zz_pX& x)

cdef extern from "NTL/lzz_pE.h" namespace "NTL":
    cdef cppclass zz_pE:
        zz_pE()

    void zz_pE_init "NTL::zz_pE::init"(const zz_pX& P)

    const zz_pX& rep(const zz_pE& a)

    cdef cppclass zz_pEPush:
        pass

cdef extern from "NTL/lzz_pEX.h" namespace "NTL":
    cdef cppclass zz_pEX:
        zz_pEX()

    long deg(const zz_pEX& a)
    void SetCoeff(zz_pEX& x, long i, long a)
    void SetCoeff(zz_pEX& x, long i)
    void SetX(zz_pEX& x)

cdef extern from "NTL/lzz_pEXFactoring.h" namespace "NTL":
    void FindRoot(zz_pE& x, const zz_pEX& f)
    zz_pE FindRoot(const zz_pEX& f)
