from ffisom import *
from rains import find_root_order

from sage.rings.infinity import Infinity
from sage.arith.all import is_prime
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.finite_rings.finite_field_constructor import GF

import sys
from sage.misc.misc import cputime

from kummer_nmod import algolist

# p n o c rains ellrains pari kummer
def benchmark_all(pbound = [3, 2**10], nbound = [3, 2**8], cbound = [1, Infinity], omax = Infinity, loops = 10, tmax = Infinity, prime = False, even = False, check = 0, fname = None, write = False, overwrite = False, verbose = True):
    if write:
        mode = 'w' if overwrite else 'a'
        f = open(fname, mode, 0)
    pmin, pmax = pbound
    nmin, nmax = nbound
    cmin, cmax = cbound
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
                print("p = {}, n = {}, (o = {}, c = {})".format(p, n, o, c))
            tloops = 0
            for l in xrange(loops):
                if (o > omax) or (o == p):
                    break
                t = cputime()
                try:
                    a, b = find_gens_rains(k_flint, k_flint, use_lucas = False)
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
                if (o != 2) or (o == p):
                    break
                t = cputime()
                try:
                    a, b = find_gens_rains(k_flint, k_flint, use_lucas = True)
                    tloops += cputime() - t
                except RuntimeError:
                    pass
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    assert(g.degree() == n)
                    assert(g == b.minpoly())
                if tloops > tmax:
                    break
            tconic = tloops / (l+1)
            tloops = 0
            for l in xrange(loops):
                t = cputime()
                try:
                    a, b = find_gens_ellrains(k_flint, k_flint)
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
                if c < cmin or c > cmax:
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
            tkummer = []
            for algo in algolist:
                if c < cmin or c > cmax:
                    break
                tloops = 0
                for l in xrange(loops):
                    t = cputime()
                    a, b = find_gens_kummer(k_flint, k_flint, n, algo)
                    tloops += cputime() - t
                    if check and (l == 0 or check > 1):
                        g = a.minpoly()
                        assert(g.degree() == n)
                        assert(g == b.minpoly())
                    if tloops > tmax:
                        break
                tkummer.append(tloops / (l+1))
            if write:
                f.write(("{} {} ({}, {}) {} {} {} {}"+" {}"*len(tkummer)+"\n").format(p, n, o, c, tcyclo, tconic, tell, tpari, *tkummer))
            else:
                sys.stdout.write(("{} {} ({}, {}) {} {} {} {}"+" {}"*len(tkummer)+"\n").format(p, n, o, c, tcyclo, tconic, tell, tpari, *tkummer))
    if write:
        f.close()

def benchmark_kummer(pbound = [3, 2**10], nbound = [3, 2**8], cbound = [1, Infinity], omax = Infinity, loops = 10, tmax = Infinity, prime = False, even = False, check = 0, fname = None, write = False, overwrite = False, verbose = True):
    if write:
        mode = 'w' if overwrite else 'a'
        f = open(fname, mode, 0)
    pmin, pmax = pbound
    nmin, nmax = nbound
    cmin, cmax = cbound
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
            if c < cmin or c > cmax:
                continue
            if verbose:
                print("p = {}, n = {}, (o = {}, c = {})".format(p, n, o, c))
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
            tkummer = []
            for algo in algolist:
                tloops = 0
                for l in xrange(loops):
                    t = cputime()
                    a, b = find_gens_kummer(k_flint, k_flint, n, algo)
                    tloops += cputime() - t
                    if check and (l == 0 or check > 1):
                        g = a.minpoly()
                        assert(g.degree() == n)
                        assert(g == b.minpoly())
                    if tloops > tmax:
                        break
                tkummer.append(tloops / (l+1))
            if write:
                f.write(("{} {} ({}, {}) {}"+ " {}"*len(algolist) + "\n").format(p, n, o, c, tpari, *tkummer))
            else:
                sys.stdout.write(("{} {} ({}, {}) {}"+ " {}"*len(algolist) + "\n").format(p, n, o, c, tpari, *tkummer))
    if write:
        f.close()
