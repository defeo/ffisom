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
`find_gens` and described in its docstring.
'''

from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.arith.all import gcd, lcm
from sage.rings.integer_ring import crt_basis
from sage.rings.finite_rings.constructor import GF
from sage.misc.cachefunc import cached_method
from sage.sets.set import Set
from sage.combinat.set_partition import SetPartitions
from sage.structure.factorization import Factorization
from sage.categories.cartesian_product import cartesian_product as CProd
from sage.misc.all import prod
from sage.misc.misc import cputime, walltime
from sage.rings.infinity import Infinity

class FalseConjecture(Exception):
    pass

def test_gens(p, n):
    '''
    Test routine for `find_gens`. Constructs two random
    extensions of F_p of degree n, then calls find_gens and
    tests that the returned elements have the same minimal polynomial
    over F_p, and that the polynomial has degree n.
    '''
    c, w = cputime(), walltime()
    k1 = GF(p**n, 'z1', modulus='random')
    k2 = GF(p**n, 'z2', modulus='random')
    print "Field creation: CPU %s, Wall %s" % (cputime(c), walltime(w))

    c, w = cputime(), walltime()
    a, b = find_gens(k1, k2)
    print "Rains' algorithm: CPU %s, Wall %s" % (cputime(c), walltime(w))

    P = a.minpoly()
    assert(P.degree() == n)
    assert(P(b) == 0)


def find_gens_list(klist, r = 0, verbose = False, omax = Infinity):
    '''
    Use Rain's method to find a generator of a subfield of degree r
    within the ambient fields in klist.

    The algorithm is done in two steps:

    1. Find a small integer m, such that the m-th roots of unity
       generate an extension of the subfields, together with some
       additional constraints detailed below. This is done by
       `find_root_order`.

    2. Find a uniquely defined Galois orbit of m-th Gaussian periods.
       Return arbitrary elements of these orbits.
       This is done by `find_unique_orbit`.

    Read the docstrings of `find_root_order` and `find_unique_orbit`
    to find out more on the algorithm.
    '''
    p = klist[0].characteristic()
    ngcd = gcd([k.degree() for k in klist])
    nlcm = lcm([k.degree() for k in klist])
    if r == 0:
        r = ngcd
    assert all(k.degree() % r == 0 for k in klist)

    # Find a suitable m, see doc below
    o, G = find_root_order(p, [ngcd, nlcm], r, verbose=verbose)
    if o > omax:
        raise RuntimeError

    # Construct extension of the fields, if needed
    #
    # Note: `find_unique_orbit` is horribly slow if k is not a
    # field. Some composita tricks would be welcome.
    return tuple(find_gen_with_data(k, r, o, G, verbose=verbose) for k in klist)

def find_gen(k, r = 0):
    return find_gens_list([k], r)[0]

def find_gens(k1, k2, r = 0, omax = Infinity, verbose = False):
    return find_gens_list([k1, k2], r, omax = omax, verbose = verbose)

def find_gen_with_data(k, r, o, G, verbose = False):
    p = k.characteristic()
    n = k.degree()
    ro = r*o
    o = ro // ro.gcd(n)

    if o > 1:
        if verbose:
            print "Using auxiliary extension of degree {}".format(o)
        P = GF(p**o, 'z').polynomial()
        from embed import computeR
        Pext = computeR(k.polynomial(), P)
        kext = GF(p**(n*o), 'zext', modulus=Pext)
        #kext = k.extension(P)

    # Compute the unique (up to Galois action) elements, descended to
    # k, if needed.
    if o == 1:
        u = find_unique_orbit(k, G)
    else:
        u = find_unique_orbit(kext, G)
        R = k.polynomial_ring()
        from embed import inverse_embed
        u = k(inverse_embed(R(u), k.polynomial(), R(P), R(kext.polynomial())))

    return u


def accept(p, n, r, l, e, s):
    '''
    This function is passed down to `sieve`. It accepts only if:

    (1)  the order t of p in ℤ/l^e is a multiple of s;
    (2)  gcd(nlcm, t / gcd(t, ngcd)) = 1.
    (3)  gcd(s, φ(m)/s) = 1.
    (4)  gcd(t/s, p) = 1.

    These three conditions imply the conditions (1)-(3) above assuming
    l is prime.
    '''
    # Degrees of the ambient fields
    ngcd, nlcm = n

    # Generic case
    if l != 2 and p != l:
        m = l**e
        phi = (l - 1) * l**(e-1)
        ord = Zmod(m)(p).multiplicative_order()
        return (ord % s.expand() == 0 and
                (phi // s.expand()).gcd(s.expand()) == 1 and
                (ord // (ord.gcd(ngcd))).gcd(nlcm) == 1 and
                (ord // s.expand()).gcd(p) == 1)

    # Special treatement for ℤ/2^x
    elif l == 2:
        return NotImplementedError

    else:
        return False

def accept_noglue(p, n, r, l, e, s):
    '''
    This function is passed down to `sieve`. It accepts only if:

    (1)  the order of p in ℤ/l^e is a multiple of s;
    (3') λ(l^e) / s is coprime to s (λ is the Carmichael function).

    These two conditions imply the conditions (1)-(3) above assuming
    l is prime.
    No gluing can be performed on accepted values if condition (2) above
    has to be preserved.
    '''
    # Generic case
    if l != 2 and p != l:
        m = l**e
        ord = (l - 1) * l**(e-1)
        return ((ord // s.expand()).gcd(s.expand()) == 1 and          # (3')
                all(Zmod(m)(p)**(ord // ell) != 1 for (ell, _) in s)) # (1)

    # Special treatement for ℤ/2^x
    elif l == 2:
        return ((e == 2 and p % 4 == 3) or
                (p != 2 and e - s[0][1] == 2 and     # (3')
                 Zmod(2**e)(p)**(s.expand() // 2)))  # (1)

    else:
        return False


def find_root_order(p, n, r = 0, accept = accept, verbose = True):
    '''
    Assuming that r divides n and p is prime,
    search for a small integer m such that:

    1. the order of <p> ⊂ ℤ/m* is equal to r⋅o;
    2. gcd(n, o) = 1;
    3. ℤ/m* = <p^o> × G for some G ⊂ ℤ/m*;

    then return o and a set of generators for G, together with their
    respective orders.

    The first condition implies that the m-th roots of unity generate
    a superfield of F_{p^r}.

    The second condition is not strictly necessary, but it makes
    the algorithm much more efficient by ensuring that the superfield
    is generated over F_{p^n} by a polynomial with coefficients in F_p.

    When r strictly divides n, the weaker condition:

    2'. gcd(n, r⋅o / gcd(n, r⋅o)) = 1; 

    can be used to work in a smaller superfield.

    The third condition is needed so that the construction of
    `find_unique_orbit` applies. In his paper, Rains' gives a
    sufficient condition for (3), namely that

    3'. gcd(r, φ(m)/r) = 1.

    It is easy to see that this condition is also necessary when ℤ/m*
    is cyclic, but there are easy counterexamples when it is not
    (e.g., take p=233, r=6, m=21).

    The integer m and the decomposition of ℤ/m are computed by the
    function `sieve` below, with acceptance criterion given by (1),
    (2) and (3').

    To conclude: bounds on m, hard to tell. When r is prime, m must be
    prime, and under GRH the best bound is m ∈ O(r^{2.4 + ε})
    [1]. Heuristically m ∈ O(r log r).  Pinch [2] and Rains give some
    tabulations.

    Note: this algorithm loops forever if p=2 and 8|r. Rains proposes
    two fixes in his paper, which we haven't implemented yet.

    [1]: D. R. Heath-Brown, Zero-free regions for Dirichlet L-functions, and
    the least prime in an arithmetic progression
    [2]: R. G. E. Pinch. Recognizing elements of finite fields.
    '''
    # For each prime power dividing r, find the smallest multiplier k.
    m = sieve(p, n, r, accept)

    if verbose:
         print "Using roots of unity of order {}".format(m)

    # Construct ℤ/m* as the product of the factors ℤ/f*
    R = Zmod(prod(f for f, _, _ in m))
    crt = map(R, crt_basis([f for f, _, _ in m]))
    G = [(sum(crt[:i])
          + R(-1 if f % 2 == 0 and f != 4 else (Zmod(f).unit_gens()[0]**s)) * crt[i]
          + sum(crt[i+1:]),
          t)
         for i, (f, s, t) in enumerate(m)]
    assert(all(g**e == 1 for (g, e) in G))

    ord = R(p).multiplicative_order()
    assert ord % r == 0

    return ord // r, G


def sieve(p, n, r = 0, accept = None):
    '''
    Given p, n and r, find the smallest integer m such that there
    exists an element of multiplicative order r in ℤ/m.

    This is equivalent to the condition: for each prime power r'|r,
    there exists a prime power m'|m such that r'|φ(m'). One immediate
    consequence is that r|φ(m).

    The output is a list of pairs of integers (m_i, r_i, o_i), such that

    - m_i is a prime power;
    - the m_i are pairwise coprime and m = ∏ m_i;
    - the r_i are pairwise coprime and r = ∏ r_i;
    - r_i · o_i = φ(m_i).

    The optional parameter `accept` can be passed in order to put
    additional constraints on m. If it is given, it must be a function
    satisfying:

    - it takes six arguments (p, n, r, l, e, s), where p, n, r, l and e
      are integers and s is the factorization of an integer
      (a `Factorization` object);
    - it returns a boolean;
    - let s and t be coprime, accept(p, n, r, l, e, s·t) returns `True` iff
      accept(p, n, r, l, e, s) and accept(p, n, r, l, e, t) also return `True`.

    Then, assuming m_i = l_i^e_i, with l_i prime, the output also
    satisfies

    - accept(p, n, r, l_i, e_i, r_i) returns `True` for any i.

    The obvious use case for `accept` is to ensure that a specific
    element of ℤ/m has order r (or divisible by r), instead of any
    element, e.g. p in the case of Rains' algorithm.

    ## Algorithm

    The algorithm starts by factoring r into prime powers r_i. For
    each r_i, it finds the smallest prime power m_i = l_i^e_i such
    that r_i|φ(m_i) and `accept(p, n, r, l_i, e_i, r_i)` returns `True`
    (if `accept` is provided).

    At this point, the lcm of all the m_i is an acceptable value for
    m, although not necessarily the smallest one. More generally, let

      r_{1,...,s} = r₁ · r₂ ··· r_s

    for some subset of the r_i (up to renumbering), then the lcm of
    m₁, ..., m_s is an acceptable value for r_{1,...,s}, although not
    necessarily the smallest one. Call m_{1,...,s} this optimal value.
    The algorithm goes on by computing the optimal values m_X for
    larger and larger subsets X ⊂ {r_i}, until r is reached.

    This is done by testing all the possible prime powers between the
    largest m_Y already computed for any proper subset Y ⊂ X, and the
    smallest lcm(m_W, m_Z) for any proper partition X = W ∪ Z with
    W ∩ Z = ∅.

    For any given r' = r_{1,...,s}, the algorithm looks for primes of
    the form k·r' + 1. If r_s = l_s^e_s is the largest prime factor of
    r' and if

      r₁ · r₂ ··· r_{s-1} | l_s - 1

    the algorithm also looks for powers of l_s. It is easy to show
    that the algorithm needs consider no other integer.

    ### Example

    Let r = 60 = 4·3·5

    The algorithm starts by computing

      m_{4} = 1·4 + 1 = 5,
      m_{3} = 2·3 + 1 = 7,
      m_{5} = 2·5 + 1 = 11.

    Then it considers all possible pairs of factors. For {4,3}, one
    possible value is 5·7 = 35. 2 divides 3-1, hence the algorithm
    looks for primes and powers of 3 comprised between 7 and 35. It
    finds that 13 is a suitable value, hence

      m_{4,3} = 1·12 + 1 = 13,

    and similarly

      m_{4,5} = 25, m_{3,5} = 31.

    Finally, it considers all three factors. It is useless to test
    prime powers below 31, because 3 and 5 must divide them. There are
    three possible partitions of {4,3,5} into two disjoints sets:

      lcm(m_{4}, m_{3,5}) = 155
      lcm(m_{3}, m_{4,5}) = 175
      lcm(m_{5}, m_{4,3}) = 143

    hence 143 is an acceptable value, and the algorithm needs to test
    no further. The first value tried is 1·60 + 1 = 61, and it turns
    out it is a prime, hence the algorithm returns

      [(61, 60)].

    ### Complexity

    The complexity is obviously polynomial in r. Heuristically, it
    should be sublinear, the dominating step being the factorization
    of r, but it is not so easy to prove it.

    If c is the number of primary factors of r, the main loop is
    executed 2^c times, which is clearly sublinear in r.

    At each iteration, all partitions of the current subset into two
    disjoint subsets must be considered, hence a very crude lower
    bound for this combinatorial step is

      ∑_{i=1}^{c} binom(c,i) (2^{i-1} - 1) < 3^c << e^{o(log r)}

    The most expensive operation of each cycle is the primality
    testing (and, eventually, the `accept` function). Heuristically,
    at each iteration O(log r) primality tests are needed, each with a
    polynomial cost in log r. As noted in `find_root_order`, the best
    bounds under GRH give O(r^{1.4 + ε}) primality tests, instead.

    Whatever the provable complexity is, this algorithm is extremely
    fast in practice, and can handle sizes which are way beyond the
    tractability of the other steps of Rains' algorithm.
    '''
    # Degrees of the the ambient fields.
    ngcd, nlcm = n

    # Actual extension degree within the ambient fields
    if r == 0:
        r = ngcd

    # If accept is not given, always accept
    if accept is None:
        accept = lambda p, n, r, l, e, s : True

    class factorization:
        '''
        This class represents the factorization of an integer,
        carrying one more piece piece of information along each
        primary factor: factors are pairs (p^e, o), with o an
        integer.
        '''
        def __init__(self, f):
            '''
            f must be a dictionary with entries of the form
                l : (e, o)
            '''
            self.factors = f

        def lcm(self, other):
            '''
            Compute the lcm of two factorizations. The auxiliary
            information is multiplied together:

              lcm ( (l^e, s), (l^d, t) ) = (l^max(e,d), s·t)
            '''
            lcm = self.factors.copy()
            for (l, (e, s)) in other.factors.iteritems():
                try:
                    E, S = lcm[l]
                    lcm[l] = (max(E, e), S*s)
                except KeyError:
                    lcm[l] = (e, s)
            return factorization(lcm)

        @cached_method
        def expand(self):
            'Return the integer represented by this factorization'
            return prod(map(lambda (l, (e, _)) : l**e, self.factors.items()))

        def __str__(self):
            return ' * '.join('%d^%d<--%d' % (l, e, s)
                              for (l, (e, s)) in self.factors.iteritems())

        def __repr__(self):
            return 'factorization(%s)' % repr(self.factors)

    # Represent the factorization of r as a set of primary factors
    fact = Set(list(r.factor()))
    # This dictionary holds the values r_X for each subset of `fact`
    optima = {}

    # Main loop, execute for each subset `S` ⊂ `fact`
    # It assumes subsets are enumerated by growing size
    for S in fact.subsets():
        # ignore the empty set
        if S.is_empty():
            continue

        # A Factorization object corresponding to S
        Sfact = Factorization(S.list())

        # find `L` the greatest prime in `S`
        L, E = max(S, key=lambda (l, e) : l)
        # the product of the remaining prime powers
        c = prod(l**e for l, e in S if l != L)
        # a boolean, True only if it is worth testing powers of `L` in
        # the sequel
        powers = (L-1) % c == 0
        # the product of all the factors of `S`
        s = c * L**E

        # For singletons, we don't have an upper bound
        if S.cardinality() == 1:
            start = 1
            end = None
        # For larger subsets, we compute the minimum and maximum value
        # to test
        else:
            parts = SetPartitions(S, 2)
            # Start from the largest m_X already computed for any
            # strict subset of `S`
            start = max(optima[T].expand() for T in S.subsets()
                        if (not T.is_empty()) and (T != S))
            # Stop at the smallest lcm of any 2-partition of `S`
            end = min((optima[T0].lcm(optima[T1]) for T0, T1 in parts),
                      key=lambda x : x.expand())

        # We only consider primes of the form
        #    k·s + 1
        # and powers of `L` of the form
        #    k·s + L^E
        # if `powers` is true.
        # We determine the starting `k`
        k = start // s or 1
        while True:
            m = k*s + 1
            # Once the upper bound, is reached, it is used with no
            # further test
            if end is not None and m >= end.expand():
                optima[S] = end
                break
            # Test primes
            elif m.is_prime() and accept(p, n, r, m, 1, Sfact):
                optima[S] = factorization({m : (1, s)})
                break
            # Test powers of `L`
            # notice the correction for L = 2 on the second line
            d = k*c + 1
            if (powers and d.is_power_of(L) and
                (L != 2 or E == 1 or k > 1) and
                accept(p, n, r, L, E + d.valuation(L), Sfact)):
                optima[S] = factorization({L : (E + d.valuation(L), s)})
                break

            k += 1

    # the last computed m_X is the optimum for r
    return [(l**e, s, (l-1)*l**(e-1) // s)
            for (l, (e, s)) in optima[fact].factors.iteritems()]


def find_unique_orbit(k, G):
    '''
    Given a superfield k of F_{p^n}, and a group G ⊂ ℤ/m*
    satisfying the conditions of `find_root_order`, i.e.

    1. <p> ⊂ ℤ/m* is of order r⋅o,
    2. gcd(n, o) = 1,
    3. ℤ/m* = <p^o> × G,

    return a Gaussian period of m-th roots of unity of k that is

    a. uniquely defined up to the action of Gal(k/F_p),
    b. stable under Gal(k/F_{p^r}),
    c. a generator of F_{p^r}.

    The Gaussian period is defined as follows. Let ζ be any m-th root
    of unity of k, define the Gaussian period

      η = ∑_{g∈G} ζ^g

    To prove (a), let e ∈ ℤ/m* be written as e = h⋅p^{o⋅i}, then

      ∑_{g∈G} (ζ^e)^g = (∑_{g∈G} ζ^{g⋅h})^{p^{o⋅i}} = η^{p^{o⋅i}},

    hence η does not depend on the choice of ζ.

    To prove (b), write p = h⋅p^{o⋅i}. Then p^r = h^r ∈ G, hence
    η^{p^r} = η by the same equation as above.

    An elementary proof of (c) eludes me. Rains proves that if m is
    squarefree, then η generates F_{p^r}. Rains' proof is based
    on the study of cyclotomic rings and their reduction at (p), a
    classical topic in algebraic number theory. The proof is not very
    hard, but I would like to find one by more elementary means, like
    for (a) and (b).
    There are examples where m is squareful, and yet η generates F_{p^n},
    though it is not true in general.
    One idea I tried, is to observe that if η is in a smaller
    subfield, then we have an equation

      ∑_{g ∈ G} ζ^g = ∑_{g ∈ xG} ζ^g

    for some coset xG where x = p^{o⋅i}. This gives a linear relation
    on 2⋅#G m-th roots of unity, something which I was trying to prove
    impossible. There is a good deal of literature on linear
    dependencies of roots of unity in ℂ, but I haven't found any on
    finite fields. An old result by Mann [1] tells that for such a
    relation to happen in ℂ, m needs to divide the product of all the
    primes up to 2⋅#G + 1 (which is unlikely to happen in our case). I
    have no idea of how to adapt this to finite fields, but it doesn't
    seem impossible to do.

    As an optimization that is not implemented here, notice that if
    ℤ/m* = <p> × G' for some other group G' ⊂ G, then there is an
    uniquely defined element

      Θ = ∑_{g ∈ G'} ζ^g

    in k, and η is the trace of θ over F_{p^r}.

    Complexity of this phase: Finding a root ζ is done in O(1) trials,
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
    return sum(zeta**prod(g**e for (g, _), e in zip(G, exps)).lift()
               for exps in CProd(map(lambda (_,x): range(x), G)))
