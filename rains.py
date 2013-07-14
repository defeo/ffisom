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
    assert(P == b.minpoly())

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

    # Construct extensions of k1 and k2, if needed
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

    The first two conditions imply that the m-th roots of unity
    generate a superfield of F_{p^n}, and that the extension can be
    defined by a polynomial with coefficients in F_p.

    The third condition is needed so that the construction of
    `find_unique_orbit` applies. In his paper, Rains' gives a
    sufficient condition for (3), namely that

    3'. gcd(n, φ(m)/n) = 1.

    It is easy to see that this condition is also necessary when ℤ/m*
    is cyclic, but there are easy counterexamples when it is not
    (e.g., take p=13, n=2, m=21, although this is a stupid example
    because m=7 would suffice).

    Remark that, when n is a prime power, m must be a prime power, so
    (3') is largely sufficient (just a little warning when n is a
    power of 2). However, did I mention that we also like to do group
    theory just for the sake of it?

    The straightforward generalization of (3') is: suppose

      ℤ/m* = C₁ × C₂ × ⋅⋅⋅ × C_r

    and suppose that <p^o> ⊂ C₁, then (3) is satisfied if and only if
    gcd(n, #C₁/n) = 1.

    The problem is that the decomposition of ℤ/m* into cyclic groups
    is not unique. And here I start getting sloppy.

    If we know the factorization of m (and we do), using an SNF we can
    compute a factorization of ℤ/m* into cyclic groups of orders

      d₁ | d₂ | ⋅⋅⋅ | λ(m)

    where λ(m) is the Carmichael function (the lcm of the orders of
    the cyclic groups associated to the prime-power factors of m). The
    SNF also gives a basis for such decomposition.

    Unless the SNF has two entries equal to λ(m) (a very rare case,
    but look at m=15⋅16, if you like nasty examples), there is a
    unique group, call it C₁, of order λ(m). Suppose that the order of
    <p> is greater than the order of any other group than C₁, then
    necessarily <p> ⊂ C₁.

    My intuition is that if these things do not happen, then a smaller
    factor m would also have worked. However, if by running this
    algorithm we bump into some nasty example, I will be even more
    pleased to look at it.

    To conclude: bounds on m, hard to tell. When n is prime, m must be
    prime, and under GRH the best bound is m ∈ O(n^{2.4 + ε})
    [1]. Heuristically m ∈ O(n log n).  Pinch [2] and Rains give some
    tabulations.

    Note: this algorithm is very naive and very slow (in some cases,
    slower than `find_unique_orbit`). Any optimizations welcome.

    [1]: D. R. Heath-Brown, Zero-free regions for Dirichlet L-functions, and
    the least prime in an arithmetic progression
    [2]: R. G. E. Pinch. Recognizing elements of finite fields.
    '''
    for m in sieve(p, n):
        if m % p == 0:
            m += 1
            continue
        R = Zmod(m)
        ord = R(p).multiplicative_order()

        # Check conditions (1) and (2)
        if ord % n == 0 and n.gcd(ord // n) == 1:

            # Get the SNF of ℤ/m*
            gens = R.unit_gens()
            A = diagonal_matrix([g.multiplicative_order() for g in gens])
            S, _, Q = A.smith_form()

            # Check that my "intuition" was correct
            if A.ncols() > 1 and ord <= S.diagonal()[-2]:
                raise FalseConjecture('Wow! Report these: p=%d, n=%d, m=%d' % (p,n,m))

            # Check condition (3).  Notice that we can divide λ(m) by
            # n⋅o, rather than n, because we know that gcd(n,o) = 1
            carmich = S.diagonal()[-1]
            if n.gcd(carmich // ord) == 1:

                # Get the new generators
                gens = [(prod(g**e for g, e in zip(gens, r)), o)
                        for r, o in zip((Q**-1).rows(), S.diagonal())]
                
                # Replace the last generator by g^n
                gens[-1] = (lambda (g, o): (g**n, o // n))(gens[-1])
                return ord // n, gens


def sieve(p, n):
    '''
    '''
    class prime_data:
        def __init__(self, prime, mul, p):
            self.prime = prime
            self.mul = mul
            self.fact = prime**mul
            self.dict = {0: (1,1)}
            self.p = p

        def __getitem__(self, k):
            if k not in self.dict:
                self.dict[k] = None

                # Primes of the form m = k⋅n + 1.  k must be prime to
                # n, in order to be able to satisfy (1), (2) and (3').
                if k % self.prime != 0:
                    m = k*fact + 1
                    if m.is_prime() and m != self.p:
                        p = Zmod(m)(self.p)
                        if p^(self.fact // self.prime) != 1:
                            ord = p.multiplicative_order()
                            self.dict[k] = (m, ord)

                # Prime powers of the form prime^a
                elif ((k == self.prime and k != 2 and self.p != self.prime) 
                      or (k == 2 and self.fact == 2 and self.p % 4 == 3)
                      or (k == 4 and self.prime == 2 and self.p != 2)):
                    m = k*fact
                    p = Zmod(m)(self.p)
                    if p^(self.fact // self.prime) != 1: 
                        ord = p.multiplicative_order()
                        self.dict[k] = (m, ord)

            return self.dict[k]
    

    facts = [[prime_data(f, m, p), 1, False]
             for (f, m) in sorted(n.factor(), key=lambda (x,y): x**y)]
    while True:
        pass
    

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

