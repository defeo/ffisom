"""
  The bench routine takes three parameters:
    p : size of the base field
    start : minimum degree (default 2)
    stop : maximum degree (default 200)

  and return a list of lists containing timings (in seconds) for
  various operations (create, embed, project, iso) for composita of
  degree n x (n+1) for n between start and stop.

  Two optional parameters 'number' and 'repeat' allow influencing the
  timeit machine (for lift and push only). See 

    sage: ?sage_timeit
"""

from sage.rings.finite_rings.constructor import GF
from sage.misc.sage_timeit import sage_timeit
from sage.rings.arith import factor
from sage.rings.ff_compositum.all import *
from sage.misc.prandom import randrange
from sage.functions.other import ceil

################# Chebyshev

def _torsion_poly(ell, P=None):
    """
    Computes the ell-th gauss period. If `P` is given, it must be a
    polynomial ring into which the result is coerced.

    This is my favourite equality:
    
    sage: all(_torsion_poly(n)(I) == I^n*lucas_number2(n,1,-1) for n in range(1,10))
    True
    """
    if P is None:
        P, R, = PolynomialRing(ZZ, 'x'), ZZ, 
    elif P.characteristic() == 0:
        R = ZZ
    else:
        R = Zp(P.characteristic(), prec=1, type='capped-rel')
    
    t = [1, 0]
    for k in range(1, ell/2 + 1):
        m = R(ell - 2*k + 2) * R(ell - 2*k + 1) / (R(ell - k) * R(k))
        t.append(-t[-2] * m)
        t.append(0)

    return P(list(reversed(t))).shift(ell % 2 - 1)


# Montgomery ladder for Pell conics
def _pellmul(x, n):
    # The zero point and ours
    A, B = 2, x
    for c in reversed(n.digits(2)):
        if c == 0:
            A, B = A**2 - 2, A*B - x
        else:
            A, B = A*B - x, B**2 - 2
    return A

def cheby(K, ell, exp):
    p = K.characteristic()
    eta = K(1)
    if p != 2:
        o = (p + 1) // ell
        while (eta**2 - 4).is_square() or _pellmul(eta, o) == 2:
            eta = K.random_element()

    P = PolynomialRing(K, 'x')
    return _torsion_poly(ell**exp, P) - eta

############ Other irreducibles

def kummer(K, ell, exp):
    eta = K.one()
    try:
        while eta.nth_root(ell):
            eta += 1
    except ValueError:
        pass
    return K.polynomial_ring().gen()**(ell**exp) - eta

def rand_irred(K, n):
    return K.polynomial_ring().irreducible_element(n)

############ Combine them together

def comp_prod(P, Q):
    n = Q.degree()
    K = P.base_ring()
    A = PolynomialRing(K, 'X')
    X = A.gen()
    AA = PolynomialRing(K, 'Y,Z')
    Y, Z = AA.gens()
    return P(Y).resultant(AA(Y**n * Q(Z/Y)), Y)(1,X)

def irred(K, n):
    F = []
    p = K.characteristic()
    for ell, exp in factor(n):
        if ell == 2:
            F.append(rand_irred(K, ell**exp))
        elif ell.divides(p-1):
            F.append(kummer(K, ell, exp))
        elif ell.divides(p+1):
            F.append(cheby(K, ell, exp))
        else:
            F.append(rand_irred(K, ell**exp))

    return reduce(comp_prod, F)

####

def test_irred(p, i1, i2, number=0, repeat=3):
    K = GF(p)
    P = FFDesc(p, irred(K, i1).list())
    Q = FFDesc(p, irred(K, i2).list())
    R = P.compositum(Q)
    
    context = globals()
    context.update(locals())
    # must be done only once to avoid (dis)counting caching
    tcomp = sage_timeit('P.compositum(Q)', context, number=1, repeat=1, seconds=True)

    a = FFElt(P, [randrange(p) for i in range(i1)])
    b = FFElt(Q, [randrange(p) for i in range(i2)])
    c = FFElt(R, [randrange(p) for i in range(i1*i2)])
    context = globals()
    context.update(locals())
    tembed_pre = sage_timeit('a.embed(b,R)', context, number=1, repeat=1, seconds=True)
    tembed = sage_timeit('a.embed(b,R)', context, number=number, repeat=repeat, seconds=True)
    
    # Burn R's precomputed polynomials
    ab = a.embed(b,R)
    ab.project(P,b)
    
    ab = a.embed(b, R)
    d = a.embed(b, R)
    context = globals()
    context.update(locals())
    tproject_pre = sage_timeit('ab.project(P,b)', context, number=1, repeat=1, seconds=True)
    tproject = sage_timeit('ab.project(P,b)', context, number=number, repeat=repeat, seconds=True)
    tmono2dual = (sage_timeit('c.embed(P, P)', context, number=1, repeat=1, seconds=True) -
                  sage_timeit('c.embed(P, P)', context, number=1, repeat=1, seconds=True))
    tmul = sage_timeit('c*c', context, number=number, repeat=repeat, seconds=True)
    ttmul = sage_timeit('c*d', context, number=number, repeat=repeat, seconds=True)
    tiso = sage_timeit('c.eltseq_dual_python(P,Q)', context, number=number, repeat=repeat, seconds=True)
    tbsgs = sage_timeit('c.eltseq_mono_BSGS(P,Q)', context, number=number, repeat=repeat, seconds=True)

    return tcomp, tembed_pre, tembed, tproject_pre, tproject, tmono2dual, tmul, ttmul, tiso, tbsgs


def bench(p, start=2, stop=200):
    l = []
    i = start
    while i < stop:
	print i,
        try:
            l.append((i, test_irred(p, i, i+1)))
        except:
            print "failed"
	else:
	    print sum(l[-1][1])
        i += ceil(i/10)
    return l

def gplot_out(l):
    s = "#\t" + "\t\t".join(('Comp\t', 'EmbedPre', 'Embed\t', 
                             'ProjectPre', 'Project\t', 'M2D\t', 
                             'Mulmod\t', 'Tmulmod\t', 'Iso\t', 'Bsgs')) + "\n"
    for i, t in l:
        s += "%d\t%s\n" % (i, "\t".join(map(str, t)))
    return s
