#-*- encoding: utf-8 

from itertools import product

v = {}
# Don't go over n = 5 !
n = 3

#for c in ['a', 'b', 'c', 'p', 'q', 'x', 'y', 'z']:
#    v[c] = [SR.var('%s%d' % (c, i)) for i in range(n)]
#PR.<X> = SR[]

R = PolynomialRing( QQ, ','.join('%s%d' % p for p in product('abpx', range(n))) )
for i, c in enumerate(['a', 'b', 'p', 'x']):
    v[c] = R.gens()[i*n:i*n+n]
PR.<X> = R[]

A, B = (sum(v[c][i]*X^i for i in range(n)) for c in ('a', 'b'))
P = X^n + sum(v['p'][i]*X^i for i in range(n))
_, iPprime, _ = xgcd( *map(lambda x: x.change_ring(R.fraction_field()), (P.derivative(), P)) )
#iPprime = iPprime.map_coefficients(lambda x : x.simplify_rational())



#####

def findmatrix(B, C):
    I = R.ideal(list(B))
    return matrix([c.lift(I) for c in C])

def mulmatrix(P, n):
    return P.reverse().sylvester_matrix(X^n)[0:n].transpose()

def mulmodmatrix(P, Mod):
    return matrix([(P*X^i % Mod).coefficients() for i in range(P.degree() + 1)]).transpose()

MAmod = mulmodmatrix(A, P)
Tmul = MAmod.transpose() * vector(v['x'])

# This is the matrix of transposed modular multiplication with respect
# to its "scalar" argument, i.e:
#   DotProd(x) : A → A*
#                a ↦ a · x
# It is a Hankel matrix
DotProd = findmatrix(v['a'], Tmul)
assert DotProd == matrix([[Tmul[i].coefficient(v['a'][j]) for j in range(n)] for i in range(n)])

####

def Tr(P):
    return P.derivative().reverse()/PowerSeriesRing(X.base_ring(), 'X', default_prec=n)(P.reverse())

def pol2tr(P, Mod):
    return mulmodmatrix(P, Mod).transpose()*vector(Tr(Mod))

def tr2pol(tr, Mod, iMod):
    return ((PR(list(tr)) * Mod.reverse()).truncate(Mod.degree()).reverse() * iMod % Mod)#.map_coefficients(lambda x : x.simplify_rational())

ATr = pol2tr(A, P)
A2 = tr2pol(ATr, P, iPprime)
assert A == A2

BTr = pol2tr(B, P)
ABTr = pol2tr(A*B % P, P)

# This is just transposed multiplication
TransMul = findmatrix(BTr, ABTr)
assert TransMul == MAmod.transpose()

########

Mpol2tr = findmatrix(v['a'], ATr)
Mtr2pol = (Mpol2tr^-1 * P.discriminant()).change_ring(R)
assert Mtr2pol.denominator() == 1

# This is the matrix of multiplication in the dual basis (up to a denominator).
# Pretty messy, huh?
DualMul = DotProd * Mtr2pol

# this verification is just too slow (blame polynomial multiplication)
assert n > 3 or (DualMul * ATr == P.discriminant() * pol2tr(tr2pol(v['x'], P, iPprime) * A % P, P))
