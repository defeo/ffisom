from sage.rings.finite_rings.constructor import GF
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.constructor import EllipticCurve

def find_gens(k1, k2, k = None, bound = None):
    '''
    Compute elements of k1 and k2 with the same minimal polynomial using
    an elliptic curve variation of Rains' algorithm.

    INPUT:

    - ``k1`` -- an extension of degree r of k.

    - ``k2`` -- an extension of degree r of k.

    - ``k`` -- (default : None) the base field.

    - ``bound`` -- (default : None) maximal value for torsion.

    EXAMPLES::

        sage: from ellrains import find_gens
        sage: R.<X> = PolynomialRing(GF(5))
        sage: f = X^19 + X^16 + 3*X^15 + 4*X^14 + 3*X^12 + 3*X^9 + 2*X^8 + 2*X^7 + 2*X^4 + X^3 + 4*X^2 + 4*X + 2
        sage: g = X^19 + 2*X^18 + 2*X^17 + 4*X^16 + X^15 + 3*X^14 + 2*X^13 + X^12 + 2*X^11 + 2*X^10 + 2*X^9 + X^8 + 4*X^6 + X^5 + 3*X^4 + 2*X^2 + 4*X + 4
        sage: k1 = GF(5**19, name = 'x', modulus = f)
        sage: k2 = GF(5**19, name = 'y', modulus = g)
        sage: tuple = find_gens(k1, k2)
        sage: tuple[0].minpoly() == tuple[1].minpoly()
        True
    '''
    if k is None:
        k = k1.base_ring()
    p = k.characteristic()
    q = k.cardinality()

    if k1.degree() != k2.degree():
        raise NotImplementedError
    r = k1.degree()
    
    # List of candidates for l.
    lT = find_l(k, r, bound)
    if lT is None:
        raise RuntimeError, "no suitable l found"

    # Find an elliptic curve with the given trace. 
    E = find_elliptic_curve(k, lT)
    if E is None:
        raise RuntimeError, "no suitable elliptic curve found"

    Ek1 = E.change_ring(k1)
    Ek2 = E.change_ring(k2)

    return (find_unique_orbit(Ek1, lT[0], r), find_unique_orbit(Ek2, lT[0], r))

def find_unique_orbit(E, l, r):
    '''
    Compute an element uniquely definied by E and l (up to Galois action).
    INPUT:
    
    - ``E`` -- an elliptic curve with a unique l-torsion subgroup.

    - ``l`` -- a prime number.

    - ``r`` -- the extension degree.

    EXAMPLES::

        sage: from ellrains import find_unique_orbit
        sage: p = 5
        sage: r = 19
        sage: l = 229
        sage: E = EllipticCurve(j = GF(p)(1))
        sage: EK = E.change_ring(GF(p**r, prefix = 'z', conway = True))
        sage: elem1 = find_unique_orbit(EK, l, r)
        sage: elem2 = find_unique_orbit(EK, l, r)
        sage: elem1.minpoly() == elem2.minpoly()
        True
    '''
    from finite_field_flint_fq_nmod import FiniteField_flint_fq_nmod
    K = E.base_ring()
    p = K.characteristic()

    if p in [2, 3]:
        raise NotImplementedError

    if r % 2 == 0:
        raise NotImplementedError

    if E.j_invariant() in [K(0), K(1728)]:
        raise NotImplementedError

    if isinstance(K, FiniteField_flint_fq_nmod):
        from xz_coordinates_flint_fq_nmod import mul_ltr
    else:
        from xz_coordinates import mul_ltr

    # Primitive root of unity of order (l-1)/r
    zeta = GF(l).multiplicative_generator()**r
    # Torsion point of order l
    P = find_torsion_point(E, mul_ltr, l)

    period = sum(xaff(mul_ltr(P, (zeta**i).lift(), E.a4(), E.a6())) for i in 
               xrange(ZZ((l-1)/(2*r))))

    return period#, P, E, l

def find_elliptic_curve(k, lT):
    '''
    Find elliptic curve with given trace mod l.

    INPUT:

    - ``k`` -- a base field,

    - ``lT`` -- a pair of an integer and a set of candidates for the trace.

    EXAMPLES::
    
        sage: from ellrains import find_elliptic_curve
        sage: lT = (229,[2])
        sage: find_elliptic_curve(GF(5), lT)
        Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 5
    '''
    l, T = lT
    L = GF(l)

    for j in k:
        # Don't mess with complicated cases
        if j == k(0) or j == k(1728):
            continue

        E = EllipticCurve(j = j)
        t = L(E.trace_of_frobenius())
        if t in T:
            return E
        elif -t in T:
            return E.quadratic_twist()

    return None

def find_traces(k, r, l, rl = None):
    '''
    Compute traces such that the charpoly of the Frobenius has a root of order n
    modulo l.

    INPUT:

    - ``k`` -- the base field.

    - ``r`` -- the degree of the extension.

    - ``l`` -- a prime number.

    EXAMPLES::

        sage: from ellrains import find_traces
        sage: r = 281
        sage: l = 3373
        sage: k = GF(1747)
        sage: find_traces(k, r, l)
        [4, 14, 18, 43, 57, 3325, 3337, 3348, 3354, 3357, 3364]
    '''
    L = GF(l)
    q = L(k.order())
    bound = 4*k.characteristic()

    if rl is None:
        rl = ZZ((l-1)/r)

    # Primitive root of unity of order r
    zeta = L.multiplicative_generator()**(ZZ((l-1)/r))
    # Candidate eigenvalue
    lmbd = L(1)
    # Traces
    T = []
    for i in xrange(r-1):
        lmbd *= zeta
        # Ensure that the order is r
        if r.gcd(i) != 1:
            continue
        trc = lmbd + q/lmbd
        # Check Hasse bound
        if trc.centerlift()**2 < bound:
           T.append(trc)

    return T

def find_l(k, r, bound = None):
    '''
    Compute a prime l such that r divides (l-1) and the corresponding traces.
 
    INPUT:

    - ``k`` -- the base field.

    - ``l`` -- the extension degree.

    - ``bound`` -- (default: None) the maximal value for (l-1)/r.

    EXAMPLES::

        sage: from ellrains import find_l
        sage: find_l(GF(1747), 281)
        (3373, [4, 14, 18, 43, 57, 3325, 3337, 3348, 3354, 3357, 3364]
        sage: find_l(GF(11), 23)
        (47, [0, 1, 2, 3, 44, 45, 46])
    '''
    q = k.order()

    if bound is None:
        bound = r**2 # Arbitrary value.  

    for rl in xrange(bound):
        # gcd(r, (l-1)/r) == 1
        if r.gcd(rl) != 1:
           continue

        l = rl*r + 1
        # l prime
        if not l.is_prime():
            continue

        # Orders of eigenvalues not dividing each other
        if r % GF(l)(q).multiplicative_order() == 0:
            continue

        T = find_traces(k, r, l)
        # At least one candidate
        if len(T) ==  0:   
            continue       
        
        return l, T

# Assume l is prime
def find_torsion_point(E, mul_ltr, l):
    K = E.base_ring()
    cofactor = E.cardinality()//l

    while True:
        x = K.random_element()
        if not E.is_x_coord(x):
            continue
        P = mul_ltr((x, K(1)), cofactor, E.a4(), E.a6())
        if P[1] == 0:
            continue
        return P

def xaff(P):
    return (P[0])/(P[1])

