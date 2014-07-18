# -*- coding: utf-8 -*-
from sage.misc.misc import cputime, walltime
from sage.rings.arith import euler_phi
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.integer_ring import ZZ
from sage.functions.other import sqrt
from sage.schemes.elliptic_curves.constructor import EllipticCurve
import XZ


def isom_elliptic(k1, k2, k = None, bound = None):
    '''
    INPUT : 
    
    - ``k1`` -- a finite field, extension of degree n of k.

    - ``k2`` -- a finite field, extension of degree n of k; k1 != k2.

    - ``k`` -- (default : None) a  finite field of characteristic p!=2, it plays
      the role of the base field for k1 & k2.

    - ``Y_coordinates`` -- (default : False) a boolean that will be used to 
      know whether to use Y or X coordinates for the Gaussian elliptic periods.

    - ``bound`` -- (default : None) a positive integer used as the max for m.

    OUTPUT : 
    
    - A tuple of unique elements with the same minimal polynomial in k1 and k2 
      respectively

    EXAMPLES :

        sage : R.<X> = PolynomialRing(GF(5))

        sage : f = X^19 + X^16 + 3*X^15 + 4*X^14 + 3*X^12 + 3*X^9 + 2*X^8 + 2*X^7 + 2*X^4 + X^3 + 4*X^2 + 4*X + 2)
        
        sage : g = X^19 + 2*X^18 + 2*X^17 + 4*X^16 + X^15 + 3*X^14 + 2*X^13 + X^12 + 2*X^11 + 2*X^10 + 2*X^9 + X^8 + 4*X^6 + X^5 + 3*X^4 + 2*X^2 + 4*X + 4

        sage : k1 = GF(5**19, name = 'x', modulus = f)

        sage : k2 = GF(5**19, name = 'y', modulus = g)

        sage : tuple = isom_elliptic(k1, k2)

        sage : tuple[0].minpoly() == tuple[1].minpoly()

        True

        sage : tuple_Y = isom_elliptic(k1, k2, Y_coordinates = True)

        sage : tuple_Y[0].minpoly() == tuple_Y[1].minpoly()

        True

    ALGORITHM:

    Given two extensions of the same base field and the same degree, we return 
    two unique elements the use of Gaussian elliptic period on normal elements,
    via the function find_unique_orbit_elliptic, on an curve E which is 
    determined by the function find_elliptic_curve.

    First we have to find an integer m and an elliptic curve E  over k such 
    that :
    TODO 

    .. TODO::

        The case j = 1728 and j = 0.
    '''
    if k is None:
	    k = k1.base_ring()
    p = k.characteristic()
    n = k1.degree()  
    q = k.cardinality()
    
    # We compute a list of candidates for m (i.e. such that n divides phi(m) 
    # and (phi(m)/n,n) = 1. It lacks the conditions on the trace.
    w = cputime()
    m_t = find_m(n, k, bound)
    w_m = cputime(w)
    
    if m_t is None:
        raise RuntimeError, 'No suitable m found, increase your bound.'

    # Finding the elliptic curve on which we can work. 
    w = cputime()
    E, case, compteur = find_elliptic_curve(k, k1, m_t) 
    w_E = cputime(w)

    if E is None:
        raise RuntimeError, 'No suitable elliptic curve found, check your \
                                                                    parameters'

    Ek1 = E.change_ring(k1)
    Ek2 = E.change_ring(k2)

    liste = (find_unique_orbit_elliptic(Ek1, m_t[0], case), 
        find_unique_orbit_elliptic(Ek2, m_t[0], case))

    w_ordm = liste[0][0]
    w_period = liste[0][1]
    return w_m, w_E, w_ordm, w_period, m_t[0]
    #return a, b

def find_unique_orbit_elliptic(E, m, case = 0):
    '''
    INPUT : 
    
    - ``E`` -- an elliptic curve with the properties given in isom_elliptic 
      and/or find_elliptic_curve.

    - ``m`` -- an integer with the properties given in isom_elliptic and/or in 
      find_m.

    - ``Y_coordinates`` -- boolean (default : False) determines if X or Y 
      coordinates are used.

    - ``case`` -- integer (default : 0) depends on the j-invariant's value :
        - ``0`` means j is not 0 nor 1728 or t = 0,
        - ``1`` means j is 1728,
        - ``2`` means j is 0.

    OUPUT : 
    
    - A uniquely defined element of the field where E is defined, namely the 
      extension of degree n considered; unique means every produced elements 
      have the same minimal polynomial.

    EXAMPLES :

    - Case j != 0, 1728

        sage: E = EllipticCurve(j = GF(5)(1))

        sage: EK = E.change_ring(GF(5**19, prefix = 'z', conway = True)

        sage: m = 229

        sage: elem1 = find_unique_orbit_elliptic(EK,m)

        sage: elem2 = find_unique_orbit_elliptic(EK,m)

        sage: elem1.minpoly() == elem2.minpoly()

        True

    - Case j = 1728 and trace != 0

        sage : E = EllipticCurve(j = GF(5)(1728))

        sage: EK = E.change_ring(GF(5**19, prefix = 'z', conway = True)

        sage: m = 229

        sage: elem1 = find_unique_orbit_elliptic(EK,m)

        sage: elem2 = find_unique_orbit_elliptic(EK,m)

        sage: elem1.minpoly() == elem2.minpoly()

        True

    - Case j = 0 and trace != 0

        sage: E = EllipticCurve(j = GF(7)(0))

        sage: EK = E.change_ring(GF(7**23, prefix = 'z', conway = True)

        sage: m = 139

        sage: elem1 = find_unique_orbit_elliptic(EK,m)

        sage: elem2 = find_unique_orbit_elliptic(EK,m)

        sage: elem1.minpoly() == elem2.minpoly()

        True


    ALGORITHM:
    TODO
    '''
    n = E.base_ring().degree()

    # Searching for a point of order exactly m.
    w = cputime()
    P = XZ.find_ordm(E, m)
    w_ordm = cputime(w)
    w = cputime()

    if case == 0:
        # Looking for a generator of order exactly phi(m)/n in 
        # phi(m)/something.
        gen_G = Integers(m).unit_gens()[0]**n
        order = euler_phi(m)//(2*n)

        r = sum((ZZ(gen_G**i)*P)[0] for i in range(order))
        w_period = cputime(w)
        return w_ordm, w_period

    elif case == 1:
        gen_G = Integers(m).unit_gens()[0]**n
        order = euler_phi(m)/(4*n)
        
        r = sum(((ZZ(gen_G**i)*P)[0])**2 for i in range(order))
        w_period = cputime(w)
        return w_ordm, w_period

    elif case == 2:
        gen_G = Integers(m).unit_gens()[0]**n
        order = euler_phi(m)/(6*n)

        r = sum(((ZZ(gen_G**i)*P)[0])**3 for i in range(order))
        w_period = cputime(w)
        return w_ordm, w_period



def find_elliptic_curve(k, K, m_t):
    '''
    INPUT: 

    - ``k`` -- a base field

    - ``K`` -- an extension of K of degree n.

    - ``m_t`` -- a list of tuple containing a integer and a set of candidates 
      for the trace.

    OUTPUT: 
    
    - An elliptic curve defined over k with the required properties.

    - An integer case, depending on the value of the j-invariant of said 
      elliptic curve.

    - An integer m statisfying the properties described in isom_elliptic which 
      we will be using to compute Gaussian periods.

    ..NOTE::

        The case j = 0 or 1728 are not implemented yet. They shall raise 
        NotImplementedError.

    EXAMPLES:
    
    sage: R.<X> = PolynomialRing(GF(5))

    sage: f = R.irreducible_element(19, algorithm = 'random')

    sage: K = GF(5**19, names = 'x', modulus = f)

    sage: m_t = [(229,{0, 1, 3})]

    sage: find_elliptic(GF(5), K, m_t)

    (Elliptic Curve defined by y^2 = x^3 + x + 2 over Finite Field of size 5,
    0,
    229)

    ALGORITHM:

    TODO : Doc is deprecated, to be redone.

    Function that finds an elliptic curve with the required charateristics, 
    those given in the function isom_elliptic.

    First, we have to determine if m is composite, a prime power or a power of 
    p, the characteristic of the base field. The first case is not implemented 
    yet. 
    We also note that the m given should satisfies several conditions based
    on the characteristic and the degree of K. See the docstrings of 
    isom_elliptic for more information.

    - If m is a power of p, the charateristic of the base field k, then we shall
      proceed as follow :

        We pick a random curve E/k and we set down t = Tr_k(Fr_E), for the
        curve to be what we want, we need :

          - t not zero,
          - (Z/m)* = <t> x S or #<t> = n 
          - m divides #E/K but not #E/L, for any intermediate extension L 
            of K/k; so we can construct points of order m such that their 
            abscissas or ordinates span exactly K. Or that we haven't any point
            of order m in any sub-extension.

        Then we test those conditions for both E and its quadratic twist, if one
        of them meet the requirements, we return it and its trace. If no 
        elliptic curves are found an error is returned.

    - If m is primer power, then we shall proceed as follow :

        We have m = l^r, for l a prime. For this method to work, we need l to 
        be an Elkies prime. A prime l is an Elkies prime for an elliptic curve 
        if the charateristic polynomial of the aforesaid elliptic curve splits 
        in GF(l).

        For now, we pick a random curve E/k and for it to work, if we set down 
        t = Tr_k(Fr_E), we need the following :

          - We have x**2 - tx + q = (x - a)(x - b) mod m, meaning the polynomial
            splits in Z/m,
          - (Z/m)* = <a> x S, with #<a> = n, 
          - ord_m(a) < ord_m(b),
          - m divides #E/K but not #E/L, for any intermediate extension L 
            of K/k; so we can construct points of order m such that their 
            abscissas or ordinates span exactly K.

        Once again, we test all that for both E and its quadratic twist; if one 
        them meet the requirements, we return it, its trace and a tuple
        containing the root a and t mod m. If none are found, there is 
        something wrong.


    - If m is composite, TODO.
    '''
    p = k.characteristic()
    q = k.cardinality()
    n = K.degree()
    m = m_t[0]
    S_t = m_t[1]
    compteur = 0

    #We start by the special cases j = 1728, 0
    E_j1728 = EllipticCurve(j = k(1728))

    if q%4 != 1:
        # If q != 1 mod 4, then there's no 4th root of unity, then magically
        # all the quartic twist are already in k and the trace is 0. We just
        # have to test the only curve y\B2 = x\B3 + x.
        compteur += 1
        if 0 in S_t:
            return E_j1728, 0, compteur
    else:
        # If q = 1 mod 4, then the trace is not 0, and we have to try four
        # trace to see which is the best candidate.
        g = k.unit_gens()[0]
        c = g**((q-1)/4)
        t = E_j1728.trace_of_frobenius()
        L = [(t*(c**i).lift(), g**i) for i in range(4)]

        for i in range(4):
            compteur += 1
            if Integers(m)(L[i][0]) in S_t:
                # E, case, t
                return E_j1728.quartic_twist(L[i][1]), 1, compteur

    E_j0 = EllipticCurve(j = k(0))

    if q%3 != 1:
        # Same as before, if q != 1 mod 6, there's no 6th root of unity in
        # GF(q) and the trace is 0 (that's pretty quick reasoning.. :D).
        # Justification will come later.
        compteur += 1
        if 0 in S_t:
            return E_j0, 0, compteur
    else:
        g = k.unit_gens()[0]
        c = g**((q-1)/6)
        t = E_j0.trace_of_frobenius()
        L = [(t*(c**i).lift(), g**i) for i in range(6)]

        for l in L:
            if Integers(m)(l[0]) in S_t:
                return E_j0.sextic_twist(l[1]), 2, compteur

    # General case
    for j in k:
        if j == 0 or j == k(1728):
            continue

        E = EllipticCurve(j = j)
        t = E.trace_of_frobenius()
        L = [(t, E), (-t, E.quadratic_twist())]

        for l in L:
            compteur +=1
            if Integers(m)(l[0]) in S_t:
                return l[1], 0, compteur

    # If no elliptic curve has been found.
    return None, -1, -1

def find_trace(n,m,k):
    '''
    INPUT : an integer n, an integer m, a base field k

    OUTPUT : a list of integer mod m or a list of a couple of integers mod m

    Algorithm :

    If m is a power of p, then we look for class modulo m with order equal to n.
    Then, we return the list of all such class.

    If m is a power of prime different from p, we look for a in (Z/m)* such 
    that :

    - ord_m(a) < ord_m(q/a) and ord_m(a) = n,

    or

    - ord_m(q/a) < ord_a and ord_m(q/a) = n.

    And we return a + q/a.

    Here a plays the role of one of the two roots of the future characteristic 
    polynomial of the Frobenius of the elliptic curve we'll use; i.e.

    X^2 - (a + q/a)*X + a*(q/a) = X^2 - t*X + q

    if we write t = a + q/a. From that, we will pick elliptic curves which have 
    one of the t's as trace of its Frobenius.
    '''
    Zm = Integers(m)
    p = k.characteristic()
    q = k.cardinality()
    sq = sqrt(float(2*q))
    q_m = Zm(q)

    # If m is a multiple of p, then we just need the trace to be of order 
    #exactly n in (Z/m)*
    if not m.is_prime_power():
        raise NotImplementedError
    elif m%p == 0:
        sol = []
        phi_m = euler_phi(m)
        alpha = phi_m/n
        g = Zm.unit_gens()[0]

        log_t = [i*alpha for i in n.coprime_integers(n)]

        for t in [g**i for i in log_t]:
            if abs(t.centerlift()) > sq:
                continue
            else:
                sol.append(t)

        return set(sol)
    # We don't want q to be of order n or dividing n, then q/a would be of order
    # n; which is unacceptable.
    elif q_m**n == 1:
        return []
    else:
        sol = []
        phi_m = euler_phi(m)
        alpha = phi_m/phi_m.gcd(n)
        g = Zm.unit_gens()[0]
        Zphi_m = Integers(phi_m)
        
        log_a = [i*alpha for i in n.coprime_integers(n)]
        a = [g**i for i in log_a]
        log_q = q_m.log(g)

        for i in range(len(log_a)):
            diff = log_q - log_a[i]
            b = g**diff
            ord_b = diff.order()

            if ord_b <= n:
                continue
            elif abs((a[i] + b).centerlift()) > sq:
                continue
            else:
                sol.append(a[i] + b)

        return set(sol)
def find_m(n, k, bound = None):
    '''
    INPUT : an integers n, a base field k, an integer bound

    OUTPUT : an integer

    Algorithm :

    Functions that given an integer n (degree of an extension) and a bound 
    returns all the candidates m, such that :

    - n|phi(m), the euler totient function,

    - (n, phi(m)/n) = 1,

    - Another one ? Maybe for q = p^d we'd want (n,d) = 1,

    We can note that if m = r^e with (e-1,n) = 1 or e = 1, then r = a*n + 1 with
    (a,n) = 1 is a suitable form for m as then phi(m) = (a*n)(an + 1)^(e-1);

    It also works in the general case if all the prime factors of m are of the 
    form a*n + 1 with (a,n) = 1. You just have to apply that to them and 
    multiply the results.
    '''
    if bound is None:
        bound_a = 100 # Arbitrary value.  
    else:
        # if m = a*n + 1 < b, then a < (b- 1)/n.
        bound_a = (bound - 1) / n 

    sol = []

    for a in range(bound_a):
        m = a*n + 1
        # m composite not implemented yet
        if not m.is_prime_power():
            continue 
        elif (euler_phi(m)/n).gcd(n) != 1:
            continue
        else:
            S_t = find_trace(n, m, k)
            if len(S_t) < 1:   # Some time in the future we'd like to have a 
                continue       # better bound than just 1.
            else:
                return m, S_t 
