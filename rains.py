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
from sage.misc.cachefunc import cached_method
from sage.sets.set import Set
from sage.combinat.set_partition import SetPartitions
from sage.structure.factorization import Factorization
from sage.combinat.cartesian_product import CartesianProduct as CProd
from sage.misc.misc import cputime, walltime, prod

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

    2. For each of k₁ and k₂, find a uniquely defined Galois orbit of
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
    if (o > 1):
        print "Using auxiliary extension, this might be slow"
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

    The integer m and the decomposition of ℤ/m are computed by the
    function `sieve` below, with acceptance criterion given by (1),
    (2) and (3').
    
    To conclude: bounds on m, hard to tell. When n is prime, m must be
    prime, and under GRH the best bound is m ∈ O(n^{2.4 + ε})
    [1]. Heuristically m ∈ O(n log n).  Pinch [2] and Rains give some
    tabulations.

    Note: this algorithm loops forever if p=2 and 8|n. Rains proposes
    two fixes in his paper, which we haven't implemented yet.

    [1]: D. R. Heath-Brown, Zero-free regions for Dirichlet L-functions, and
    the least prime in an arithmetic progression
    [2]: R. G. E. Pinch. Recognizing elements of finite fields.
    '''
    def accept(r, e, n):
        '''
        This function is passed down to `sieve`. It accepts only if:

        (1)  the order of p in ℤ/r^e is a multiple of n;
        (3') λ(r^e) / n is coprime to n (λ is the Carmichael function).

        These two conditions are equivalent to the conditions (1)-(3)
        above when r is prime.
        '''
        # Generic case
        if r != 2 and p != r:
            m = r**e
            ord = (r - 1) * r**(e-1)
            return ((ord // n.expand()).gcd(n.expand()) == 1 and          # (3')
                    all(Zmod(m)(p)**(ord // ell) != 1 for (ell,_) in n))  # (1)

        # Special treatement for ℤ/2^x
        elif r == 2:
            return ((e == 2 and p % 4 == 3) or
                    (p != 2 and e - n[0][1] == 2 and     # (3')
                     Zmod(2**e)(p)**(n.expand() // 2)))  # (1)

        else:
            return False
    
    # For each prime power, find the smallest multiplier k.
    m = sieve(n, accept)

    # Construct ℤ/m* as the product of the factors ℤ/f*
    R = Zmod(prod(r for r, _, _ in m))
    crt = map(R, crt_basis([r for r, _, _ in m]))
    G = [(sum(crt[:i])
          + R(-1 if r % 2 == 0 and r != 4 else (Zmod(r).unit_gens()[0]**o)) * crt[i]
          + sum(crt[i+1:]), 
          c)
         for i, (r, o, c) in enumerate(m)]
    assert(all(g**e == 1 for (g,e) in G))

    ord = R(p).multiplicative_order()
        
    return ord // n, G


def sieve(n, accept=None):
    '''
    Given n, find the smallest integer m such that there exist an
    element of multiplicative order n in ℤ/m.

    This is equivalent to the condition: for each prime power q|n,
    there exists a prime prime power r|m such that q|φ(r). One
    immediate consequence is that n|φ(m).

    The output is a list of pairs of integers (r_i, o_i, c_i), such that

    - r_i is a prime power;
    - the r_i are pairwise coprime and m = ∏ r_i;
    - the o_i are pairwise coprime and n = ∏ o_i;
    - o_i · c_i = φ(r_i).

    The optional parameter `accept` can be passed in order to put
    additional constraints on m. If it is given, it must be a function
    satisfying:
    
    - it takes three arguments (p, e, o), where p and e are integers
      and o is the factorization of an integer (a `Factorization`
      object);
    - it returns a boolean;
    - let s and t be coprime, accept(p,e,s·t) returns `True` iff
      accept(p,e,s) and accept(p,e,t) also return `True`.

    Then, assuming r_i = p_i^e_i, with p_i prime, the output also
    satisfies

    - accept(p_i, e_i, o_i) returns `True` for any i.

    The obvious use case for `accept` is to insure that a spefic
    element of ℤ/m has order n (or divisible by n), instead of any
    element.


    ## Algorithm
    
    The algorithm starts by factoring n into prime powers q_i. For
    each q_i, it finds the smallest prime power r_i = p_i^e_i such
    that q_i|r_i and `accept(p_i, e_i, q_i)` returns `True` (if
    `accept` is provided).

    At this point, the lcm of all the r_i is an acceptable value for
    m, although not necessarily the smallest one. More generally, let

      n_{1,...,s} = q₁ · q₂ ··· q_s

    for some subset of the q_i (up to renumbmering), then the lcm of
    r₁, ..., r_s is an acceptable value for n_{1,...,s}, although not
    necessarily the smallest one. Call r_{1,...,s} this optimal value.
    The algorithm goes on by computing the optimal values r_X for
    larger and larger subsets X ⊂ {q_i}, until n is reached. 

    This is done by testing all the possible prime powers between the
    largest r_Y already computed for any proper subset Y ⊂ X, and the
    smallest lcm(r_W, r_Z) for any proper partition X = W ∪ Z with
    W ∩ Z = ∅.

    For any given n' = n_{1,...,s}, the algorithm looks for primes of
    the form k·n' + 1. If q_s = p_s^e_s is the largest prime factor of
    n' and if
    
      q₁ · q₂ ··· q_{s-1} | p_s - 1

    the algorithm also looks for powers of p_s. It is easy to show
    that the algorithm needs consider no other integer.

    ### Example

    Let n = 60 = 4·3·5

    The algorithm starts by computing

      r_{4} = 1·4 + 1 = 5,
      r_{3} = 2·3 + 1 = 7,
      r_{5} = 2·5 + 1 = 11.
    
    Then it considers all possible pairs of factors. For {4,3}, one
    possible value is 5·7 = 35. 2 divides 3-1, hence the algorithm
    looks for primes and powers of 3 comprised between 7 and 35. It
    finds that 13 is a suitable value, hence

      r_{4,3} = 1·12 + 1 = 13,

    and similarly

      r_{4,5} = 25, r_{3,5} = 31.
      
    Finally, it considers all three factors. It is useless to test
    prime powers below 31, because 3 and 5 must divide them. There are
    three possible partitions of {4,3,5} into two disjoints sets:

      lcm(r_{4}, r_{3,5}) = 155
      lcm(r_{3}, r_{4,5}) = 175
      lcm(r_{5}, r_{4,3}) = 143

    hence 143 is an acceptable value, and the algorithm needs to test
    no further. The first value tried is 1·60 + 1 = 61, and it turns
    out it is a prime, hence the algorithm returns

      [(61, 60)].

    ### Complexity

    The complexity is obviously polynomial in n. Heuristically, it
    should be sublinear, the dominating step being the factorization
    of n, but it is not so easy to prove it.

    If c is the number of primary factors of n, the main loop is
    executed 2^c times, which is clearly sublinear in n.

    At each iteration, all partitions of the current subset into two
    disjoint subsets must be considered, hence a very crude lower
    bound for this combinatorial step is

      ∑_{i=1}^{c} binom(c,i) (2^{i-1} - 1) < 3^c << e^{o(log n)}

    The most expensive operation of each cycle is the primality
    testing (and, eventually, the `accept` function). Heuristically,
    at each iteration O(log n) primality tests are needed, each with a
    polynomial cost in log n. As noted in `find_root_order`, the best
    bounds under GRH give O(n^{1.4 + ε}) primality tests, instead.

    Whatever the provable complexity is, this algorithm is extremely
    fast in practice, and can handle sizes which are way beyond the
    tractability of the other steps of Rains' algorithm.
    '''
    # If accept is not given, always accept
    if accept is None:
        accept = lambda p, e, o : True

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
                p : (e, o)
            '''
            self.factors = f

        def lcm(self, other):
            '''
            Compute the lcm of two factorizations. The auxiliary
            information is multiplied together:

              lcm ( (p^e, s), (p^d, t) ) = (p^max(e,d), s·t)
            '''
            lcm = self.factors.copy()
            for (p,(e,o)) in other.factors.iteritems():
                try:
                    E, O = lcm[p]
                    lcm[p] = (max(E,e), O*o)
                except KeyError:
                    lcm[p] = (e, o)
            return factorization(lcm)

        @cached_method
        def expand(self):
            'Return the integer represented by this factorization'
            return prod(map(lambda (p,(e,_)) : p**e, self.factors.items()))

        def __str__(self):
            return ' * '.join('%d^%d<--%d' % (p, e, o) 
                              for (p,(e,o)) in self.factors.iteritems())

        def __repr__(self):
            return 'factorization(%s)' % repr(self.factors)
            

    # Represent the factorization of n as a set of primary factors
    fact = Set(list(n.factor()))
    # This dictionary holds the values r_X for each subset of `fact`
    optima = {}

    # Main loop, execute for each subset `S` ⊂ `fact`
    for S in fact.subsets():
        # ignore the empty set
        if S.is_empty():
            continue

        # A Factorization object corresponding to S
        Sfact = Factorization(S.list())
        
        # find `P` the greatest prime in `S`
        P, E = max(S, key=lambda (p,e) : p)
        # the product of the remaining prime powers
        c = prod(p**e for p, e in S if p != P)
        # a boolean, True only if it is worth testing powers of `P` in
        # the sequel
        powers = (P-1) % c == 0
        # the product of all the factors of `S`
        o = c * P**E
        
        # For singletons, we don't have an upper bound
        if S.cardinality() == 1:
            start = 1
            end = None
        # For larger subsets, we compute the minimum and maximum value
        # to test
        else:
            parts = SetPartitions(S, 2)
            # Start from the largest r_X already computed for any
            # subset of `S`
            start = max(optima[s].expand() for p in parts for s in p)
            # Stop at the smallest lcm of any partition of `S`
            end = min((optima[p[0]].lcm(optima[p[1]]) for p in parts),
                      key=lambda x : x.expand())
            
        # We only consider primes of the form
        #    k·o + 1
        # and powers of `P` of the form
        #    k·o + P^E
        # if `powers` is true.

        # We determine the starting `k`
        k = start // o or 1
        while True:
            m = k*o + 1
            # Once the upper bound, is reached, it is used with no
            # further test
            if end is not None and m >= end.expand():
                optima[S] = end
                break
            # Test primes
            elif m.is_prime() and accept(m, 1, Sfact):
                optima[S] = factorization({m:(1, o)})
                break
            # Test powers of `P`
            # notice the correction for P = 2 on the second line
            d = k*c + 1
            if (powers and d.is_power_of(P) and
                (P != 2 or E == 1 or k > 1) and
                accept(P, E + d.valuation(P), Sfact)):
                optima[S] = factorization({P:(E + d.valuation(P), o)})
                break
                
            k += 1
            
    # the last computed r_X is the optimum for n
    return [(p**e, o, (p-1)*p**(e-1) // o)
            for (p,(e,o)) in optima[fact].factors.iteritems()]


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
    return sum(zeta**prod(g**e for (g, _), e in zip(G, exps)).lift()
               for exps in CProd(*map(lambda (_,x): xrange(x),
                                      G)))

