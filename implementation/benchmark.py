from ffisom import *
from rains import find_root_order

from sage.rings.infinity import Infinity
from sage.arith.all import is_prime
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.finite_rings.finite_field_constructor import GF

import sys
from sage.misc.misc import cputime

# p n o c rains ellrains pari kummer
def benchmark(pbound = [5, 2**10], nbound = [5, 2**8], loops = 10, omax = Infinity, cmax = Infinity, tmax = Infinity, prime = False, even = False, check = 0, fname = None, write = False, overwrite = False, verbose = True):
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
                a, b = find_gens_kummer(k_flint, k_flint)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummer = tloops / (l+1)
            if write:
                f.write("{} {} ({}, {}) {} {} {} {}\n".format(p, n, o, c, tcyclo, tell, tpari, tkummer))
            else:
                sys.stdout.write("{} {} ({}, {}) {} {} {} {}\n".format(p, n, o, c, tcyclo, tell, tpari, tkummer))
    if write:
        f.close()

def benchmark_kummer(pbound = [5, 2**10], nbound = [5, 2**8], loops = 10, omax = Infinity, cmax = Infinity, tmax = Infinity, prime = False, even = False, check = 0, fname = None, write = False, overwrite = False, verbose = True):
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
            c = Mod(p,n).multiplicative_order()
            if c > cmax:
                continue
            if verbose:
                print p, n, c
            tloops = 0
            for l in xrange(loops):
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
                t = cputime()
                # MC, MP
                a, b = find_gens_kummer(k_flint, k_flint, n, 0, 0)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummermcmp = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                t = cputime()
                # MC, no MP
                a, b = find_gens_kummer(k_flint, k_flint, n, 0, 1<<30)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummermc = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                t = cputime()
                # LA, MP
                a, b = find_gens_kummer(k_flint, k_flint, n, 1<<30, 0)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummerlamp = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                t = cputime()
                # LA, no MP
                a, b = find_gens_kummer(k_flint, k_flint, n, 1<<30, 1<<30)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummerla = tloops / (l+1)
            if write:
                f.write("{} {} ({}, {}) {} {} {} {} {}\n".format(p, n, o, c, tpari, tkummermcmp, tkummermc, tkummerlamp, tkummerla))
            else:
                sys.stdout.write("{} {} ({}, {}) {} {} {} {} {}\n".format(p, n, o, c, tpari, tkummermcmp, tkummermc, tkummerlamp, tkummerla))
    if write:
        f.close()

def benchmark_trivial(pbound = [5, 2**10], nbound = [5, 2**10], loops = 10, tmax = Infinity, check = 0, fname = None, write = False, overwrite = False, verbose = True):
    if write:
        mode = 'w' if overwrite else 'a'
        f = open(fname, mode, 0)
    pmin, pmax = pbound
    nmin, nmax = nbound
    for p in xrange(pmin, pmax):
        p = ZZ(p)
        if not p.is_prime():
            continue
        for n in ZZ(p-1).prime_divisors():
            if n < nmin or n > nmax:
                continue
            k = GF(p**n, name='z')
            k_flint = GF_flint(p, k.modulus(), name='z')
            tloops = 0
            for l in xrange(loops):
                t = cputime()
                # MC, MP
                a, b = find_gens_kummer(k_flint, k_flint, n, 0, 0)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummermcmp = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                t = cputime()
                # LA, MP
                a, b = find_gens_kummer(k_flint, k_flint, n, 1<<30, 0)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummerla = tloops / (l+1)
            if write:
                f.write("{} {}: {} {}\n".format(p, n, tkummermcmp, tkummerla))
            else:
                sys.stdout.write("{} {}: {} {}\n".format(p, n, tkummermcmp, tkummerla))
    if write:
        f.close()

def benchmark_linalg_nontriv(pbound = [5, 2**10], nbound = [5, 2**8], loops = 10, cmax = Infinity, tmax = Infinity, prime = False, even = False, check = 0, fname = None, write = False, overwrite = False, verbose = True):
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
            c = Mod(p,n).multiplicative_order()
            if c == 1 or c > cmax:
                continue
            if verbose:
                print p, n, c
                sys.stdout.flush()
            tloops = 0
            for l in xrange(loops):
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
                t = cputime()
                # MC, MP
                a, b = find_gens_kummer(k_flint, k_flint, n, 0, 0)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummermcmp = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                t = cputime()
                # MC, no MP
                a, b = find_gens_kummer(k_flint, k_flint, n, 0, 1<<30)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummermc = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                t = cputime()
                # LA, MP
                a, b = find_gens_kummer(k_flint, k_flint, n, 1<<30, 1<<30)
                tloops += cputime() - t
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tkummerla = tloops / (l+1)
            if write:
                f.write("{} {} ({}) {} {} {} {}\n".format(p, n, c, tpari, tkummermcmp, tkummermc, tkummerla))
            else:
                sys.stdout.write("{} {} ({}) {} {} {} {}\n".format(p, n, c, tpari, tkummermcmp, tkummermc, tkummerla))
    if write:
        f.close()
