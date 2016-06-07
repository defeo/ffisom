from ffisom import *
from rains import find_root_order

from sage.rings.infinity import Infinity
from sage.arith.all import is_prime
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.finite_rings.finite_field_constructor import GF

import sys
from sage.misc.misc import cputime

# p n o rains ellrains pari javad
def benchmark(pbound = [5, 2**10], nbound = [5, 2**8], loops = 10, omax = Infinity, cmax = Infinity, tmax = Infinity, prime = False, even = False, check = 0, fname = None, write = True, overwrite = False, verbose = True):
    if write:
        mode = 'w' if overwrite else 'a'
        f = open(fname, mode, 0)
    pmin, pmax = pbound
    nmin, nmax = nbound
    for p in xrange(pmin, pmax):
        p = ZZ(p)
        if not p.is_prime():
            continue
        for n in xrange(nmin, nmax):
            n = ZZ(n)
            if prime and not is_prime(n):
                continue
            if n % p == 0:
                continue
            if (not even) and (n % 2 == 0):
                continue
            k = GF(p**n, name='z')
            k_flint = GF_flint(p, k.modulus(), name='z')
            o, _ = find_root_order(p, [n, n], n, verbose=False)
            c = Mod(p,n).multiplicative_order()
            if verbose:
                print p, n, o, c
            tloops = 0
            for l in xrange(loops):
                if (o > omax) or (o == p):
                    break
                t = cputime()
                try:
                    a, b = find_gens_cyclo(k_flint, k_flint)
                    tloops += cputime() - t
                except RuntimeError:
                    pass
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tcyclo = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                t = cputime()
                try:
                    a, b = find_gens_ell(k_flint, k_flint)
                    tloops += cputime() - t
                except RuntimeError:
                    break
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tell = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                if c > cmax:
                    break
                t = cputime()
                a, b = find_gens_pari(k, k)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tpari = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                if c > cmax:
                    break
                t = cputime()
                a, b = find_gens_javad(k_flint, k_flint)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tjavad = tloops / (l+1)
            if write:
                f.write("{} {} ({}, {}) {} {} {} {}\n".format(p, n, o, c, tcyclo, tell, tpari, tjavad))
            else:
                sys.stdout.write("{} {} ({}, {}) {} {} {} {}\n".format(p, n, o, c, tcyclo, tell, tpari, tjavad))
    if write:
        f.close()
