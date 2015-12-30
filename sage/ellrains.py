from sage.rings.integer_ring import ZZ
from sage.rings.arith import gcd, lcm
from sage.rings.finite_rings.constructor import GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve

def find_gens_list(klist, r = 0, bound = None, verbose = True):
    """
    Use an elliptic curve variation of Rain's method to find a generator
    of a subfield of degree r within the ambient fields in klist.
    """
    p = klist[0].characteristic()
    ngcd = gcd([k.degree() for k in klist])
    nlcm = lcm([k.degree() for k in klist])
    if r == 0:
        r = ngcd
    assert all(k.degree() % r == 0 for k in klist)

    # This might be useful if elliptic curves defined over an 
    # extension of F_p are used.
    kb = k.base_ring()

    # List of candidates for l.
    lT = find_l(kb, r, bound)
    if lT is None:
        raise RuntimeError, "no suitable l found"

    # Find an elliptic curve with the given trace. 
    E = find_elliptic_curve(kb, lT)
    if E is None:
        raise RuntimeError, "no suitable elliptic curve found"

    return tuple(find_unique_orbit(E.change_ring(k), lT[0], r) for k in klist)

def find_gen(k, r = 0, bound = None):
    return find_gens_list([k], r, bound)

def find_gens(k1, k2, r = 0, bound = None):
    return find_gens_list([k1, k2], r, bound)


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

