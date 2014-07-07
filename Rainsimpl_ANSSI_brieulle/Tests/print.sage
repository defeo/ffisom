from rainselliptictest import *
import sys, getopt

# Evoulution of certain parameters for a single characteristic when the degree 
# is varying.
def routine_p(p, borne_nbn, start_n = 7, leap = 1): #1
    n = start_n
    print '#n m t_E ord_m period'

    for i in range(borne_nbn):
        R = PolynomialRing(GF(p), name = 'X')
        f = R.irreducible_element(n, algorithm='random')
        g = R.irreducible_element(n, algorithm='random')
        if g == f:
            g = R.irreducible_element(n, algorithm='random') 
        k1 = GF(p**n, name = 'x', modulus = f)
        k2 = GF(p**n, name = 'y', modulus = g)

        # Times to compute m, times to find an E, times to compute a point of 
        # order m, time to compute the periods, (value of m, won't be used in 
        # this routine).
        w_m, w_E, w_ordm, w_period, m = isom_elliptic(k1, k2)
        
        # Skipping first line because of an anomaly that gives an execution time
        # higher than it should be for one of the parameters. The culprit is 
        # probably a module called the first time.
        if i:
            print '%s %s %s %s %s' % (n, w_m, w_E, w_ordm, w_period)
        
        for j in range(leap):
            n = next_prime(n)

# Evolution of other parameters with p fixed.
def routine_mp(p, borne_nbn, start_n = 3, leap = 1): #2
    n = start_n
    print '#n m len(S_t) t_E compteur'

    for i in range(borne_nbn):
        m, S_t = find_m(n, GF(p))
        k1 = GF(p**n, name = 'x', modulus =
                PolynomialRing(GF(p), name='X').irreducible_element(n,
                    algorithm='random'))
        w = cputime()
        n_E = find_elliptic_curve(GF(p), k1, (m, S_t))[-1]
        w_t = cputime(w)
        if i :
            #The degree, the value of m, the number of trace candidates,
            # time to find an E, number of elliptic curves tested.
            print '%s %s %s %s %s' % (n, m, len(S_t), w_t, n_E)

        for j in range(leap):
            n = next_prime(n)

# Evolution of parameters with n fixed.
def routine_n(n, borne_nbp, start_p = 5, leap = 1): #3
    p = start_p
    print '#p m t_E ord_m period'

    for i in range(borne_nbp):
        R = PolynomialRing(GF(p), name = 'X')
        f = R.irreducible_element(n, algorithm='random')
        g = R.irreducible_element(n, algorithm='random')
        if g == f:
            g = R.irreducible_element(n, algorithm='random')

        k1 = GF(p**n, name = 'x', modulus = f)
        k2 = GF(p**n, name = 'y', modulus = g)

        w_m, w_E, w_ordm, w_period, m = isom_elliptic(k1, k2)

        if i:
            print '%s %s %s %s %s' % (p, w_m, w_E, w_ordm, w_period)
        
        for j in range(leap):
            p = next_prime(p)

# Evolution of parameters with n fixed.
def routine_mn(n, borne_nbp, start_p = 5, leap = 1): #4
    p = start_p
    print '#p m len(S_t) t_E n_E'

    for i in range(borne_nbp):
        m, S_t = find_m(n, GF(p))
        k1 = GF(p**n, name = 'x', modulus =
                PolynomialRing(GF(p), name='X').irreducible_element(n,
                    algorithm='random'))
        w = cputime()
        n_E = find_elliptic_curve(GF(p), k1, (m, S_t))[-1]
        w_E = cputime(c)
        if i :
            print '%s %s %s %s %s' % (p, m, len(S_t), w_E, n_E)

        for j in range(leap):
            p = next_prime(p)

test = int(sys.argv[1])
if test == 1:
    routine_p(ZZ(sys.argv[2]), int(sys.argv[3]))
elif test == 2:
    routine_mp(ZZ(sys.argv[2]), int(sys.argv[3]))
elif test == 3:
    routine_n(ZZ(sys.argv[2]), int(sys.argv[3]))
elif test == 4:
    routine_mn(ZZ(sys.argv[2]), int(sys.argv[3]))
else:
    raise RuntimeError, 'There\'s only four tests for now.'
