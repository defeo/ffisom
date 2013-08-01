# -*- coding: utf-8 -*-

'''
This is an implementation of E. Rains'
<http://www.math.caltech.edu/people/rains.html> algorithm for
isomorphisms of finite fields.

Given two finite fields represented as F_p[X]/f(X) and F_p[X]/g(X),
with f and g irreducible of degree n, it finds generators α, β such
that the extension of α ↦ β defines an isomorphism, and outputs
change-of-basis matrices for the two directions.

The algorithm for finding the generators α and β is implemented by
`find_gens_cyclotomic` and described in its docstring.
'''

from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.integer_ring import crt_basis
from sage.rings.finite_rings.constructor import GF
from sage.matrix.constructor import diagonal_matrix
from sage.combinat.cartesian_product import CartesianProduct as CProd
from sage.misc.misc import cputime, walltime

class FalseConjecture(Exception):
    pass

def test_gens_cyclotomic(p, n):
    '''
    Test routine for `find_gens_cyclotomic`. Constructs two random
    extensions of F_p of degree n, then calls find_gens_cyclotomic and
    tests that the returned elements have the same minimal polynomial
    over F_p, and that the polynomial has degree n.
    '''
    c, w = cputime(), walltime()
    k1 = GF(p**n, 'z1', modulus='random')
    k2 = GF(p**n, 'z2', modulus='random')
    print "Field creation: CPU %s, Wall %s" % (cputime(c), walltime(w))

    c, w = cputime(), walltime()
    a, b = find_gens_cyclotomic(k1, k2)
    print "Rains' algorithm: CPU %s, Wall %s" % (cputime(c), walltime(w))

    P = a.minpoly()
    assert(P.degree() == n)
    assert(P(b) == 0)

def find_gens_cyclotomic(k1, k2):
    '''
    Use Rain's method to find generators α ∈ k₁ and β ∈ k₂ such that
    α ↦ β defines a valid isomorphism.

    The algorithm is done in two steps:

    1. Find a small integer m, such that the m-th roots of unity
       generate an extension of k₁ and k₂, together with some
       additional constraints detailed below. This is done by
       `find_root_order`.

    2. For each of k₁ and k₂, find an uniquely defined Galois orbit of
       m-th Gaussian periods. Return arbitrary elements of those
       orbits. This is done by `find_unique_orbit`.

    Read the docstrings of `find_root_order` and `find_unique_orbit`
    to find out more on the algorithm.
    '''
    p = k1.characteristic()
    n = k1.degree()
    assert (p, n) == (k2.characteristic(), k2.degree())

    # Find a suitable m, see doc below
    o, G = find_root_order(p, n)

    # Construct extensions of k₁ and k₂, if needed
    #
    # Note: `find_unique_orbit` is horribly slow if k₁ and k₂ are not
    # fields. Some composita tricks would be welcome.
    P = None
    if (ord > n):
        P = GF(p**o, 'z').polynomial()
        k1, k2 = [k.extension(P) for k in (k1, k2)]
        
    # Return the unique (up to Galois action) elements, descended to
    # k₁ and k₂, if needed.
    return map(lambda x: x if P is None else x[0],
               [find_unique_orbit(k, G) for k in (k1, k2)])



def find_root_order(p, n):
    '''
    Search for a small integer m such that:

    1. the order of <p> ⊂ ℤ/m* is equal to n⋅o;
    2. gcd(n, o) = 1;
    3. ℤ/m* = <p^o> × G for some G ⊂ ℤ/m*;

    then return o and a set of generators for G, together with their
    respective orders.

    Rains also requires m to be square-free, but we like to live
    dangerously.

    The first condition implies that the m-th roots of unity generate
    a superfield of F_{p^n}. The second condition is not strictly
    necessary, but it makes the algorithm much more efficient by
    insuring that the superfield is generated over F_{p^n} by a
    polynomial with coefficients in F_p.

    The third condition is needed so that the construction of
    `find_unique_orbit` applies. In his paper, Rains' gives a
    sufficient condition for (3), namely that

    3'. gcd(n, φ(m)/n) = 1.

    It is easy to see that this condition is also necessary when ℤ/m*
    is cyclic, but there are easy counterexamples when it is not
    (e.g., take p=233, n=6, m=21).

    The algorithm implemented here factors n, then for each prime
    power finds an m_i satisfying (1)-(3). Since the m_i must be prime
    powers, we use condition (3') (with some adaptation for powers of
    2) to insure (3). Then, the lcm of all the m_i is used as m.
    
    This algorithm always computes the smallest m when n is a prime
    power. It may compute a subotimal m otherwise: e.g., p=2, n=6
    gives m=21, but m=9 is also a solution. It is doable, although
    painful to adapt the algorithm to always find the smallest m.

    To conclude: bounds on m, hard to tell. When n is prime, m must be
    prime, and under GRH the best bound is m ∈ O(n^{2.4 + ε})
    [1]. Heuristically m ∈ O(n log n).  Pinch [2] and Rains give some
    tabulations.

    [1]: D. R. Heath-Brown, Zero-free regions for Dirichlet L-functions, and
    the least prime in an arithmetic progression
    [2]: R. G. E. Pinch. Recognizing elements of finite fields.
    '''
    def accept(k, ell, e, p):
        '''
        This function tells wether there is a m ~ k·ℓ^e such that p
        and ℓ^e satisfy conditions (1)-(3).

        It returns m as a pair of a prime and an exponent
        '''
        n = ell**e

        # Case 1: m prime
        #
        # Primes of the form m = k⋅n + 1.  k must be prime to
        # n, in order to be able to satisfy (1), (2) and (3').
        if k % ell != 0:
            m = k*n + 1
            if m.is_prime() and m != p:
                p = Zmod(m)(p)
                ord = m - 1
                if p**(ord // ell) != 1:
                    return m, 1
                    
        # Case 2: m a higher power of a prime
        #
        # These cases are avoided by the original algorithm (which
        # needs m squarefree).
        #
        # Notice that in order to satisfy (1)-(3), there is
        # essentially only one choice for m.

        # m = ℓ^{e+1}, ℓ odd.
        elif (k == ell and k != 2 and p != ell):
            m = k*n
            p = Zmod(m)(p)
            ord = (ell - 1) * n
            if p**(ord // ell) != 1:
                return ell, e+1

        # m = 4. This only works for n=2, and hence 
        #   p = -1 mod 4
        elif (k == 2 and n == 2 and p % 4 == 3):
            return ell, 2

        # m = 2^{e+2}.
        elif (k == 4 and ell == 2 and p != 2):
            m = 4*n
            R = Zmod(m)
            p = R(p)
            if p**(n // 2) != 1:
                return ell, e+2

        return None
    
    class factor:
        '''
        This class represents one of the factors m_i. It contains
        three fields:
        
        - A prime `ell`;
        - An exponent `mul`, m_i is defined as ell^mul;
        - An integer `exp` dividing φ(m_i), it is factor of n.
        '''
        def __init__(self, ell):
            self.ell = ell
            self.mul = 1
            self.exp = 1

        def add(self, mul, exp):
            '''
            This function is called when two factors of n collide (the
            m_i's are powers of the same prime). This function updates
            self with the new data.
            '''
            self.mul = max(self.mul, mul)
            self.exp *= exp
            return self
            
        def fact(self):
            '''
            Returns m_i
            '''
            return self.ell**self.mul

        def order(self):
            '''
            Returns φ(m_i).
            '''
            return (self.ell-1) * self.ell**(self.mul-1)

        def gen(self):
            '''
            Returns a generator of the unique subgroup of order
            φ(m_i)/exp of ℤ/m_i.
            '''
            R = Zmod(self.fact())
            if R.order() == 4:
                return R(1)
            elif R.order() % 2 == 0:
                return R(-1)
            else:
                return R.unit_gens()[0]**self.exp

    # For each prime power, find the smallest multiplier k.
    sieve = {}
    for ell, mul in n.factor():
        k, m = 0, None
        while m is None:
            k += 1
            m = accept(k, ell, mul, p)
        sieve[m[0]] = sieve.get(m[0], factor(m[0])).add(m[1], ell**mul)

    # Construct ℤ/m* as the product of the factors ℤ/f*
    m = prod(sieve[ell].fact() for ell in sieve)
    R = Zmod(m)
    crt = crt_basis([sieve[ell].fact() for ell in sieve])
    G = [(sum(crt[:i]) + R(sieve[ell].gen())*crt[i] + sum(crt[i+1:]), 
          sieve[ell].order() // sieve[ell].exp)
         for (i, ell) in enumerate(sieve)]
    assert(all(g**e == 1 for (g,e) in G))

    ord = R(p).multiplicative_order()
        
    return ord // n, G
    


def find_unique_orbit(k, G):
    '''
    Given a field k isomorphic to F_{p^{n⋅o}}, and a group G ⊂ ℤ/m*
    satisfying the conditions of `find_root_order`, i.e.

    1. <p> ⊂ ℤ/m* is of order n⋅o,
    2. gcd(n, o) = 1,
    3. ℤ/m* = <p^o> × G,

    return a Gaussian period of m-th roots of unity of k that is

    a. uniquely defined up to the action of Gal(k/F_p),
    b. stable under Gal(k/F_{p^n}),
    c. a generator of F_{p^n}.

    The Gaussian period is defined as follows. Let ζ be any m-th root
    of unity of k, define the Gaussian period

      η = ∑_{g∈G} ζ^g

    To prove (a), let n ∈ ℤ/m* be written as n = h⋅p^{o⋅i}, then

      ∑_{g∈G} (ζ^n)^g = (∑_{g∈G} ζ^{g⋅h})^{p^{o⋅i}} = η^{p^{o⋅i}}.

    To prove (b), write p = h⋅p^{o⋅i}. Then p^n = h^n ∈ G, hence
    η^{p^n} = η by the same equation as above.

    An elementary proof of (c) eludes me. Rains proves that if m is
    squarefree, then η generates F_{p^n}. There are examples where m
    is squareful, and yet η generates F_{p^n}. Rains' proof is based
    on the study of cyclotomic rings and their reduction at (p), a
    classical topic in algebraic number theory. The proof is not very
    hard, but I would like to find one by more elementary means, like
    for (a) and (b).

    One idea I tried, is to observe that if η is in a smaller
    subfield, then we have an equation

      ∑_{g ∈ G} ζ^g = ∑_{g ∈ xG} ζ^g

    for some coset xG where x = p^{o⋅i}. This gives a linear relation
    on 2⋅#G m-th roots of unity, something which I was trying to prove
    impossible. There is a good deal of literature on linear
    dependencies of roots of unity in ℂ, but I haven't find any on
    finite fields. An old result by Mann [1] tells that for such a
    relation to happen in ℂ, m needs to divide the product of all the
    primes up to 2⋅#G + 1 (which is unlikely to happen in our case). I
    have no idea of how to adapt this to finite fields, but it doesn't
    seem impossible to do.

    As an optimization that is not implemented here, notice that if
    ℤ/m* = <p> × G' for some other group G' ⊂ G, then there is an
    uniquely defined element

      Θ = ∑_{g ∈ G'} ζ^g

    in k, and η is the trace of θ over F_{p^n}.

    Complexity of this phase. Finding a root ζ is done in O(1) trials,
    each costing O(n⋅o) multiplications in k. The Gaussian period is
    the sum of O(m/n) elements of k, each of which is computed with
    O(log m) multiplications in k. By computing the powers ζ^g smartly
    (not done here), one can obviously do all the job with only
    O(log m) multiplications and O(m/n) additions.

    [1]: H.B. Mann. On linear relations between roots of unity. 1965
    '''
    m = G[0][0].parent().order()
    cofactor = k.cardinality() // m
    fact = m.factor()

    # find an m-th root of unity
    zeta = 1
    while any(zeta**(m // f[0]) == 1 for f in fact):
        zeta = k.random_element()**cofactor
        
    # return the Gaussian period
    # ... lovely combinatorial iterators!
    return sum(zeta**prod(g**e for (g, _), e in zip(G, exps))
               for exps in CProd(*map(lambda (_,x): xrange(x),
                                      G)))

