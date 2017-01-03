from sage.rings.infinity import Infinity
from sage.arith.all import is_prime, is_prime_power
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.finite_rings.finite_field_constructor import GF

import os, errno, sys, time
# cputime might get screwed up by openblas
from sage.misc.misc import walltime, cputime
mytime = walltime
from cysignals.alarm import alarm, AlarmInterrupt, cancel_alarm

from ffisom import *
from sage.interfaces.magma import Magma

# p n (o, c) magma rains conic elliptic pari kummer
def benchmark(pbound = [3, 2**10], nbound = [3, 2**8], cbound = [1, Infinity], obound = [1, Infinity], loops = 10, tloop = Infinity, tmax = Infinity, prime = False, even = False, check = False, fname = None, write = False, overwrite = False, verbose = True, skip_pari = False, skip_magma = False, skip_rains = False, skip_kummer = False):
    if write:
        mode = 'w' if overwrite else 'a'
        f = open(fname, mode, 0)
    else:
        f = sys.stdout
    pmin, pmax = pbound
    nmin, nmax = nbound
    omin, omax = obound
    cmin, cmax = cbound
    M = Magma()
    for p in xrange(pmin, pmax):
        p = ZZ(p)
        if not p.is_prime():
            continue
        for n in xrange(nmin, nmax):
            n = ZZ(n)
            if (prime == 1 and not is_prime(n)) or (prime == 2 and not is_prime_power(n)):
                continue
            if n < 2:
                continue
            if n % p == 0:
                continue
            if (not even) and (n % 2 == 0):
                continue
            o, G = find_root_order(p, [n, n], n, verbose=False)
            m = G[0][0].parent().order()
            c = Mod(p,n).multiplicative_order()
            if verbose:
                sys.stdout.write("\r"+" "*79)
                print("\rp = {}, n = {}, (o = {}, c = {})".format(p, n, o, c))
            if verbose:
                t = mytime()
                sys.stdout.write("Constructing fields ({})".format(time.strftime("%c")))
                sys.stdout.flush()
            q = p**n
            k = GF(q, name='z')
            k_rand = GF(q, modulus='random', name='z')
            k_flint = GF_flint(p, k.modulus(), name='z')
            if verbose > 1:
                sys.stdout.write("\ntotal: {}s\n".format(mytime(t)))
                sys.stdout.flush()
            # Magma
            if verbose:
                sys.stdout.write("\r"+" "*79)
                sys.stdout.write("\rMagma ({})".format(time.strftime("%c")))
                sys.stdout.flush()
            tloops = 0
            for l in xrange(loops):
                if skip_magma:
                    break
                if (o > omax) or (o == p):
                    break
                # let's assume that launching a new Magma instance is cheaper
                # than computing random irreducible polynomials
                try:
                    M._start()
                except OSError as err:
                    # but it can also cause fork issues...
                    # let's accept this
                    # and fail as the situation will only worsen
                    # unless it is "just" a memory issue
                    # which should be mitigated by COW but is not
                    #print(err)
                    if err.errno == errno.ENOMEM:
                        break
                    else:
                        raise
                try:
                    k_magma = M(k)
                    k_rand_magma = M(k_rand)
                    if tloop is not Infinity:
                        alarm(tloop)
                    t = mytime()
                    k_magma.Embed(k_rand_magma, nvals=0)
                    #M._eval_line("Embed(k_magma, k_rand_magma);", wait_for_prompt=False)
                    tloops += mytime(t)
                except TypeError:
                    # sage/magma interface sometimes gets confused
                    pass
                except (KeyboardInterrupt, AlarmInterrupt):
                    # sage interface eats KeyboardInterrupt
                    # and AlarmInterrupt derives from it
                    tloops = 0
                    break
                finally:
                    if tloop is not Infinity:
                        cancel_alarm()
                    M.quit()
                    # sage pexpect interface leaves zombies around
                    try:
                        while os.waitpid(-1, os.WNOHANG)[0]:
                            pass
                    # but sometimes every child is already buried
                    # and we get an ECHILD error...
                    except OSError:
                        pass
                if tloops > tmax:
                    break
            tmagma = tloops / (l+1)
            if verbose > 1:
                sys.stdout.write("\ntotal: {}s, per loop: {}s\n".format(tloops, tloops/(l+1)))
                sys.stdout.flush()
            # Rains algorithms
            if verbose:
                sys.stdout.write("\r"+" "*79)
                sys.stdout.write("\rCyclotomic Rains ({})".format(time.strftime("%c")))
                sys.stdout.flush()
            trains = []
            tloops = 0
            for l in xrange(loops):
                if skip_rains:
                    break
                if (o > omax) or (o == p):
                    break
                try:
                    if tloop is not Infinity:
                        alarm(tloop)
                    t = mytime()
                    a, b = find_gens_cyclorains(k_flint, k_flint, use_lucas = False)
                    tloops += mytime(t)
                except AlarmInterrupt:
                    tloops = 0
                    break
                finally:
                    if tloop is not Infinity:
                        cancel_alarm()
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    if g.degree() != n:
                        raise RuntimeError("wrong degree")
                    if g != b.minpoly():
                        raise RuntimeError("different minpolys")
                if tloops > tmax:
                    break
            trains.append(tloops / (l+1))
            if verbose > 1:
                sys.stdout.write("\ntotal: {}s, per loop: {}s\n".format(tloops, tloops/(l+1)))
                sys.stdout.flush()
            # Conic Rains
            if verbose:
                sys.stdout.write("\r"+" "*79)
                sys.stdout.write("\rConic Rains ({})".format(time.strftime("%c")))
                sys.stdout.flush()
            tloops = 0
            for l in xrange(loops):
                if skip_rains:
                    break
                if (o != 2) or (o > omax) or (o == p):
                    break
                try:
                    if tloop is not Infinity:
                        alarm(tloop)
                    t = mytime()
                    a, b = find_gens_cyclorains(k_flint, k_flint, use_lucas = True)
                    tloops += mytime(t)
                except AlarmInterrupt:
                    tloops = 0
                    break
                finally:
                    if tloop is not Infinity:
                        cancel_alarm()
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    if g.degree() != n:
                        raise RuntimeError("wrong degree")
                    if g != b.minpoly():
                        raise RuntimeError("different minpolys")
                if tloops > tmax:
                    break
            trains.append(tloops / (l+1))
            if verbose > 1:
                sys.stdout.write("\ntotal: {}s, per loop: {}s\n".format(tloops, tloops/(l+1)))
                sys.stdout.flush()
            # Elliptic Rains
            if verbose:
                sys.stdout.write("\r"+" "*79)
                sys.stdout.write("\rElliptic Rains ({})".format(time.strftime("%c")))
                sys.stdout.flush()
            tloops = 0
            for l in xrange(loops):
                if skip_rains:
                    break
                try:
                    if tloop is not Infinity:
                        alarm(tloop)
                    t = mytime()
                    a, b = find_gens_ellrains(k_flint, k_flint)
                    tloops += mytime(t)
                except RuntimeError:
                    # sometimes no suitable elliptic curve exists
                    pass
                except AlarmInterrupt:
                    tloops = 0
                    break
                finally:
                    if tloop is not Infinity:
                        cancel_alarm()
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    if g.degree() != n:
                        raise RuntimeError("wrong degree")
                    if g != b.minpoly():
                        raise RuntimeError("different minpolys")
                if tloops > tmax:
                    break
            trains.append(tloops / (l+1))
            if verbose > 1:
                sys.stdout.write("\ntotal: {}s, per loop: {}s\n".format(tloops, tloops/(l+1)))
                sys.stdout.flush()
            # PARI/GP
            if verbose:
                sys.stdout.write("\r"+" "*79)
                sys.stdout.write("\rPARI/GP ({})".format(time.strftime("%c")))
                sys.stdout.flush()
            tloops = 0
            for l in xrange(loops):
                if skip_pari:
                    break
                if c < cmin or c > cmax:
                    break
                try:
                    if tloop is not Infinity:
                        alarm(tloop)
                    t = mytime()
                    a, b = find_gens_pari(k, k)
                    tloops += mytime(t)
                except AlarmInterrupt:
                    tloops = 0
                    break
                finally:
                    if tloop is not Infinity:
                        cancel_alarm()
                if check and (l == 0 or check > 1):
                    g = a.minpoly()
                    if g.degree() != n:
                        raise RuntimeError("wrong degree")
                    if g != b.minpoly():
                        raise RuntimeError("different minpolys")
                if tloops > tmax:
                    break
            tpari = tloops / (l+1)
            # Kummer algorithms
            tkummer = []
            # only linalg and modcomp implemented for c==1
            for i, algo in enumerate(kummer_algolist[2*(c==1):-2*(c==1)-1]):
                if verbose:
                    sys.stdout.write("\r"+" "*79)
                    sys.stdout.write("\rKummer {} ({})".format(kummer_namelist[2*(c==1)+i], time.strftime("%c")))
                    sys.stdout.flush()
                tloops = 0
                for l in xrange(loops):
                    if skip_kummer:
                        break
                    if c < cmin or c > cmax:
                        break
                    try:
                        if tloop is not Infinity:
                            alarm(tloop)
                        t = mytime()
                        a, b = find_gens_kummer(k_flint, k_flint, n, algo)
                        tloops += mytime(t)
                    except AlarmInterrupt:
                        tloops = 0
                        break
                    finally:
                        if tloop is not Infinity:
                            cancel_alarm()
                    if check and (l == 0 or check > 1):
                        g = a.minpoly()
                        if g.degree() != n:
                            raise RuntimeError("wrong degree")
                        if g != b.minpoly():
                            raise RuntimeError("different minpolys")
                    if tloops > tmax:
                        break
                tkummer.append(tloops / (l+1))
                if verbose > 1:
                    sys.stdout.write("\ntotal: {}s, per loop: {}s\n".format(tloops, tloops/(l+1)))
                    sys.stdout.flush()
            tkummer = 2*(c == 1)*[0] + tkummer + 2*(c == 1)*[0]
            if verbose:
                sys.stdout.write("\r")
                sys.stdout.flush()
            f.write(("{} {} ({}, {})" + " {}" + " {}"*len(trains) + " {}" + " {}"*len(tkummer)+"\n").format(p, n, o, c, *([tmagma] + trains + [tpari] + tkummer)))
    if write:
        f.close()
