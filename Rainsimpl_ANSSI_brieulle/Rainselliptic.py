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

    ALGORITHM:

    Given two extensions of the same degree defined by two different
    polynomials, we want to find, with the help of elliptic period, two elements
    with an unique orbit under the action of the Galois group. The algorithm is
    as follows :

        - First we have to find an integer m with the following properties :

            - We want to have n | phi(m) and (phi(m)/n, n) = 1.
            - We need m such that there exist an eigenvalue of the Frobenius of
              order n in (Z/m)* and for that egeinvalue to be of minimal order. 
              From there we can construct a good class for the trace of the 
              Frobenius, if we have one or more of those classes, we select 
              this m (note : for now, m is just a prime or a prime power).

          This is done by the function find_m.
    
        - After that, we need to pick a good elliptic curve. Which is an 
          elliptic curve defined over GF(q) with its trace of Frobenius has its 
          class modulo m equal to one of the good candidates, the one returned 
          by the function find_m. The properties on m ensure that the elliptic 
          curve over GF(q^n) has points of order m and that the abscissas of 
          such points span GF(q^n). The function doing that is 
          find_elliptic_curve.

        - Finally, we compute the elliptic periods u1 and u2 on both E/k1 and 
          E/k2 using the abscissas of a point of order m on each curves, 
          the isomorphism we are looking for is the one sending u1 on u2. 
    '''
    if k is None:
	    k = k1.base_ring()
    p = k.characteristic()
    n = k1.degree()  
    q = k.cardinality()
    
    # We compute a list of candidates for m (i.e. such that n divides phi(m) 
    # and (phi(m)/n,n) = 1. It lacks the conditions on the trace.
    m_t = find_m(n, k, bound)
    
    if m_t is None:
        raise RuntimeError, 'No suitable m found, increase your bound.'

    # Finding the elliptic curve on which we can work. 
    E, case = find_elliptic_curve(k, k1, m_t) 

    if E is None:
        raise RuntimeError, 'No suitable elliptic curve found, check your \
                                                                    parameters'

    Ek1 = E.change_ring(k1)
    Ek2 = E.change_ring(k2)

    a, b = (find_unique_orbit_elliptic(Ek1, m_t[0], 
        case), find_unique_orbit_elliptic(Ek2, m_t[0], case))

    return a, b

def find_unique_orbit_elliptic(E, m, case = 0):
    '''
    INPUT : 
    
    - ``E`` -- an elliptic curve with the properties given in isom_elliptic 
      and/or find_elliptic_curve.

    - ``m`` -- an integer with the properties given in isom_elliptic and/or in 
      find_m.

    - ``case`` -- integer (default : 0) depends on the j-invariant's value :
        - ``0`` means j is not 0 nor 1728 or E is supersingular,
        - ``1`` means j is 1728,
        - ``2`` means j is 0.

    OUPUT : 
    
    - An element in the field K_E over which E is defined, with a unique orbit 
      under the action of the Galois group  K_E/k.

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

    From a point of order m on E/GF(q^n), we use its abscissas to generate a 
    uniquely defined element. To defined such element, we need to calculate 
    periods of the Galois action. The trace of the elliptic curve we are using 
    is of the form t = a + q/a, with a of order n in (Z/m)*. So for S a subgroup
    of the Galois groupe, we have (Z/m)* = <a> x S. To compute the elliptic
    periods, we use the formulas :

        - u = sum_{i \in S} (([i]P)[0])^2, for j not 0 nor 1728 or t = 0,
        - u = sum_{i \in S} (([i]P)[0])^4, for j = 1728,
        - u = sum_{i \in S} (([i]P)[0])^6, for j = 0.
    '''
    n = E.base_ring().degree()

    # Loking for a point of order exactly m.
    P = XZ.find_ordm(E, m)

    if case == 0:
        # Looking for a generator of order exactly phi(m)/n in 
        # (Z/m)*/something.
        gen_G = Integers(m).unit_gens()[0]**n
        order = euler_phi(m)//(2*n)

        return sum((ZZ(gen_G**i)*P)[0] for i in range(order))
    elif case == 1:
        gen_G = Integers(m).unit_gens()[0]**n
        order = euler_phi(m)/(4*n)
        
        return sum(((ZZ(gen_G**i)*P)[0])**2 for i in range(order))

    elif case == 2:
        gen_G = Integers(m).unit_gens()[0]**n
        order = euler_phi(m)/(6*n)

        return sum(((ZZ(gen_G**i)*P)[0])**3 for i in range(order))

def find_elliptic_curve(k, K, m_t):
    '''
    INPUT: 

    - ``k`` -- a base field,

    - ``K`` -- an extension of K of degree n,

    - ``m_t`` -- a list of tuple containing an integer and a set of candidates 
      for the trace.

    OUTPUT: 
    
    - An elliptic curve defined over k with the required properties.

    - An integer case, depending on the value of the j-invariant or the 
      the supersingularity of said elliptic curve.

    EXAMPLES:
    
    sage: R.<X> = PolynomialRing(GF(5))

    sage: f = R.irreducible_element(19, algorithm = 'random')

    sage: K = GF(5**19, names = 'x', modulus = f)

    sage: m_t = (229,{2})

    sage: find_elliptic_curve(GF(5), K, m_t)

    (Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 5, 1)

    ALGORITHM:

    The goal is to pick an elliptic curve of which the trace of its Frobenius 
    t is among the right class modulo m, the ones in S_t. We do that in order to 
    have point of order m only on E/GF(q^n) or above, since then the abscissas 
    of a point of order m will span GF(q^n) and we'll be able to compute the 
    elliptic periods as it was planned.

    We start by looking at the two special cases j = 0 and 1728. If q is not 1 
    modulo 4 and 3 respectively, then the curves are supersingular (t = 0) and 
    if 0 is among the good traces, they are to be treated like the other 
    curves.

    If for j = 0 we have q = 1 mod 3, then we have to tests E(j = 0) and all of
    its sextic twists. Once again if t is in S_t, then we return the right 
    curves and the case 2 to compute the periods accordingly.

    If for j = 1728 we have q = 1 mod 4, then we have to tests  E(j = 1728) and
    all of its quartic twists. If one the trace is in S_t, we return the right 
    curves and the case 1.

    If j != 0 and 1728, then we tests all the elements of GF(q) to find the 
    right curve. For each j we test if t or -t is in S_t and return the curve 
    accordingly plus the case 0.

    If no curves are found, we return None and the case -1, which will raise an
    runtimeError in the main code.
    '''
    p = k.characteristic()
    q = k.cardinality()
    n = K.degree()
    m = m_t[0]
    S_t = m_t[1]

    #We start by the special cases j = 1728, 0
    E_j1728 = EllipticCurve(j = k(1728))

    if q%4 != 1:
        # If q != 1 mod 4, then there's no 4th root of unity, then magically
        # all the quartic twist are already in k and the trace is 0. We just
        # have to test the only curve y² = x³ + x.
        if 0 in S_t:
            return E_j1728, 0
    else:
        # If q = 1 mod 4, then the trace is not 0, and we have to try four
        # trace to see which is the best candidate.
        g = k.unit_gens()[0]
        c = g**((q-1)/4)
        t = E_j1728.trace_of_frobenius()
        L = [(t*(c**i).lift(), g**i) for i in range(4)]

        for i in range(4):
            if Integers(m)(L[i][0]) in S_t:
                # E, case, t
                return E_j1728.quartic_twist(L[i][1]), 1

    E_j0 = EllipticCurve(j = k(0))

    if q%3 != 1:
        # Same as before, if q != 1 mod 6, there's no 6th root of unity in
        # GF(q) and the trace is 0 (that's pretty quick reasoning.. :D).
        # Justification will come later. Since q = 1 mod 2, if q = 1 mod 3
        # then q = 1 mod 6.
        if 0 in S_t:
            return E_j0, 0
    else:
        g = k.unit_gens()[0]
        c = g**((q-1)/6)
        t = E_j0.trace_of_frobenius()
        L = [(t*(c**i).lift(), g**i) for i in range(6)]

        for l in L:
            if Integers(m)(l[0]) in S_t:
                return E_j0.sextic_twist(l[1]), 2

    # General case
    for j in k:
        if j == 0 or j == k(1728):
            continue

        E = EllipticCurve(j = j)
        t = E.trace_of_frobenius()
        L = [(t, E), (-t, E.quadratic_twist())]

        for l in L:
            if Integers(m)(l[0]) in S_t:
                return l[1], 0

    # If no elliptic curve has been found.
    return None, -1

def find_trace(n,m,k):
    '''
    INPUT :

    - ``n`` -- an integer, the degree of the extension,

    - ``m`` -- an integer, a candidate for the paramater m,

    - ``k`` -- a finite field, the base field.

    OUTPUT : 

    - A list of integer modulo m with the good properties.

    EXAMPLES :

    sage: n = 281

    sage: m = 3373

    sage: k = GF(1747)

    sage: find_trace(n,m,k)

    {4, 14, 18, 43, 57, 3325, 3337, 3348, 3354, 3357, 3364}

    ALGORITHM :

    The algorithm is pretty straightforward. We select all the elements of 
    order n and look for some properties of their class modulo m and add 
    them to list if they meet the requirements. They will be the class 
    candidates for the trace of the future elliptic curves.

    - If m is a power of p, the characteristic, then we look for all the 
    elements of order n and check if they end in the Hasse interval. 

    - If m is a prime (power) different from p, then we start by computing the 
    logarithm of elements of order n in (Z/m)*. You have the minimal polynomial
    of the Frobenius equal to X**2 - t*X + q = (X - a)(X - q/a) mod m. We look 
    for a among the element of order n modulo m such that the other root q/a is
    of greater order. If their sum a + q/a = t falls into the Hasse interval, 
    then we add t in the good candidates.

    - If m is composite, we raise a NotImplementedError.
    '''
    Zm = Integers(m)
    p = k.characteristic()
    q = k.cardinality()
    sq = sqrt(float(2*q))
    q_m = Zm(q)

    # This test may be obsolete if the function is used inside 'isom_elliptic'.
    if not m.is_prime_power():
        raise NotImplementedError
    # If m is a power of p, then we just need the trace to be of order 
    # exactly n in (Z/m)*
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
        
        # Computing the logarithm of element of order n.
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
    INPUT :

    - ``n`` -- an integer, the degree,

    - ``k`` -- a finite field, the base field,

    - ``bound`` -- (default : None) a positive integer used as the max for m.

    OUTPUT :

    - A tuple containing an integer and a set of class modulo m.

    EXAMPLES :

    sage: find_m(281, GF(1747))

    (3373, {4, 14, 18, 43, 57, 3325, 3337, 3348, 3354, 3357, 3364})

    sage: find_m(23, GF(11))

    (47, {0, 1, 2, 3, 44, 45, 46})

    ALGORITHM :
    
    First and foremost we are looking for an integer m for which n | phi(m). A 
    good way to obtain such integers is to look for those of the form 

    m = an + 1, 

    because then phi(m) = d.(an) which is divisible by n. We also want phi(m) 
    to be coprime with n, then choosing m to be a prime (which is possible 
    thanks to Dirichlet's theorem on the arithmetic progressions) we ensure 
    that it is actually the case.

    It still works fine, theoratically, if an + 1 is a prime power. Though, we 
    almost get to pick a m that is prime.
    
    Once we have that integer, we try to compute good candidates for the 
    traces and see how many works. If less than a certain number works (this 
    number is equal to 1 at the moment), we discard it and test the next prime 
    power. When one is found, we return it with its trace class candidates.
    '''
    if bound is None:
        bound_a = 100 # Arbitrary value.  
    else:
        # if m = a*n + 1 < b, then a < (b- 1)/n.
        bound_a = (bound - 1) / n 

    for a in range(bound_a):
        m = a*n + 1
        # m composite not implemented yet
        if not m.is_prime_power():
            continue 
        elif (euler_phi(m)/n).gcd(n) != 1:
            continue
        else:
            S_t = find_trace(n, m, k)
            # Some time in the future we'd like to have a 
            # better bound than just 1.
            if len(S_t) < 1:   
                continue       
            else:
                return m, S_t 
