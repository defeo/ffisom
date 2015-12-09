from sage.misc.sage_timeit import sage_timeit
from sage.rings.arith import inverse_mod

############################################################
############################################################
##
## some functions for finite field embeddings and the like
## linear forms are represented as polynomials
##
############################################################
############################################################



############################################################
############################################################
## Basic univariate polynomial operations
############################################################
############################################################

############################################################
## 1/rev(P) mod x^m
############################################################
def S(P, m):
    PP = P.parent()
    from sage.rings.laurent_series_ring import LaurentSeriesRing
    LS = LaurentSeriesRing(PP.base_ring(), name='z', default_prec=m)
    irP = 1/LS(P.reverse())
    return PP(irP.list())

############################################################
## particular case of the previous one, needed for tmulmod
############################################################
def S0(P):
    return S(P, P.degree()-1)

############################################################
## computes the k first Newton sums of M
## (as a polynomial)
############################################################
def trace_from_minpoly(M, k):
    m = M.degree()
    tmp = M.derivative().reverse(m-1)*S(M, k) ## m-1 = degree
    return tmp.truncate(k)

############################################################
## composed product of P and Q
############################################################
def computeR(P, Q):
    m = P.degree()
    n = Q.degree()
    t = trace_from_minpoly(P, 2*m*n+4)
    u = trace_from_minpoly(Q, 2*m*n+4)
    v = [t[i]*u[i] for i in range(2*m*n+4)]
    from sage.matrix.berlekamp_massey import berlekamp_massey
    BM = berlekamp_massey(v)
    return BM


############################################################
############################################################
## Transposed algorithms
############################################################
############################################################

############################################################
## univariate transposed product
## K^{m+k} -> K^k, with deg(B) <= m
############################################################
def tmul(A, B, m, k):
    return ((A*(B.reverse(m))).truncate(k+m)).shift(-m) ## m = degree

############################################################
## transposed remainder. ell is reduced mod P, outputs k terms
############################################################
def trem(ell, P, k):
    m = P.degree()
    U = P.parent()
    G = S(P, k-m)
    A = tmul(ell, P, m, k-m)
    C = (G*A).truncate(k-m)
    return ell - C.shift(m)

############################################################
## transposed modular multiplication, B.ell mod P
## S = S0(P), precomputed
############################################################
def tmulmod(ell, B, P, S):
    m = P.degree()
    U = P.parent()
    A = tmul(U(ell), P, m, m-1)
    C = (S*A).truncate(m-1)
    D = U(ell) - C.shift(m)
    return tmul(D, B, m-1, m)


############################################################
## transpose of the right multiplication of a polynomial matrix Lin
## by a polynomial matrix B. After transposition, the
## elementary products are K^{m+k} -> K^k, with deg(B[i,ell]) <= m
############################################################
def trightmatmul(Lin, B, m, k):
    L = Lin.transpose()
    qq, rr = B.dimensions()
    rr, pp = L.dimensions()
    C = matrix(B[0,0].parent(), qq, pp)
    for i in range(qq):
        for j in range(pp):
            for ell in range(rr):
                C[i,j] = C[i,j] + tmul(L[ell,j], B[i,ell], m, k)
    return C.transpose()


############################################################
############################################################
## Trace formulas
############################################################
############################################################

############################################################
## given traces tr(A B^i) and the min poly M of B,
## recovers the expression A = C(B)
############################################################
def convert_from_trace(t, M):
    m = M.degree()
    D = inverse_mod(M.derivative(), M)
    N = (M.reverse(m)*t).truncate(m)
    Nstar = N.reverse(m-1)
    C = (Nstar*D) % M
    return C

############################################################
## bivariate conversion from trace
## t is a bivariate polynomial
############################################################
def convert_from_trace_bi(t, P, Q):
    M = t.parent()
    _x, _y = M.gens()[:1]
    m = P.degree()
    n = Q.degree()
    IP = inverse_mod(P.derivative(), P)
    IQ = inverse_mod(Q.derivative(), Q)
    N = (P.reverse(m)(_x) * Q.reverse(n)(_y) * t).truncate(_x, m).truncate(_y, n)
    Nstar = N(1/_x,1/_y) * _x^(m-1) * _y^(n-1)
    Nstar = Nstar.numerator()
    H = (Nstar*IP(_x)*IQ(_y)).reduce([P(_x), Q(_y)])
    return H


############################################################
############################################################
## Embedding and its inverse
############################################################
############################################################


############################################################
## TODO: for most functions, we should use precomputations for
## many quantities (traces, S0, ...)
############################################################


############################################################
## embeds F(X) mod R
## SP = S0(P), precomputed
############################################################
def embed(F, P, Q, R):
    U = P.parent()
    m = P.degree()
    n = Q.degree()

    SP = S0(P)
    t = trace_from_minpoly(P, m)
    u = trace_from_minpoly(Q, m*n)
    ell = tmulmod(t, F, P, SP)
    ellstar = trem(ell, P, m*n)
    v = U([ellstar[i]*u[i] for i in range(m*n)])
    return convert_from_trace(v, R)

############################################################
## transpose of the previous one
############################################################
def tembed(ell, P, Q, R):
    U = P.parent()
    m = P.degree()
    n = Q.degree()
    SP = S0(P)
    SR = S0(R)
    t = trace_from_minpoly(P, m)
    u = trace_from_minpoly(Q, m*n)
    v = convert_from_trace(ell, R)
    ellstar = U([v[i]*u[i] for i in range(m*n)])
    ellS = ellstar % P
    return tmulmod(t, ellS, P, SP)

############################################################
## inverse embedding
############################################################
def inverse_embed(G, P, Q, R):
    m = P.degree()
    n = Q.degree()
    SP = S0(P)
    SR = S0(R)
    t = trace_from_minpoly(R, m*n)
    ell = tmulmod(t, G, R, SR)
    v = tembed(ell, P, Q, R)
    return convert_from_trace(v, P)/n

############################################################
############################################################
## Change of basis and its inverse
############################################################
############################################################

############################################################
## change of basis, 1st algorithm
############################################################
def change_basis(FF, P, Q, R):
    M = FF.parent()
    U = P.parent()
    m = P.degree()

    SP = S0(P)
    SR = S0(R)
    SQ = S0(Q)

    S = embed(U.gen(), P, Q, R)
    G = 0
    for i in range(m-1,-1,-1):
        Fi = U(FF.coefficient({M.gens()[0]:i})(0,U.gen()))
        Gi = embed(Fi, Q, P, R)
        G = (G*S+Gi) % R
    return G

############################################################
## transposed change of basis, 1st algorithm
## M is the target bivariate polynomial ring
############################################################
def tchange_basis(GG, P, Q, R, M):
    U = P.parent()
    m = P.degree()

    SP = S0(P)
    SR = S0(R)
    SQ = S0(Q)

    S = embed(U.gen(), P, Q, R)
    ell = M(0)
    for i in range(m):
        tmp = tembed(GG, Q, P, R)
        ell = ell + tmp(M.gens()[1])*M.gens()[0]^i
        GG = tmulmod(GG, S, R, SR)
    return ell

############################################################
## change of basis, second algorithm
############################################################
def change_basis_2(FF, P, Q, R):
    M = FF.parent()
    U = P.parent()
    m = P.degree()
    n = Q.degree()

    SP = S0(P)
    SR = S0(R)
    SQ = S0(Q)

    np = n + m - 1
    p = ceil(sqrt(np))
    q = ceil(np/p)
    T = embed(U.gen(), Q, P, R)
    iT = inverse_mod(T, R)
    iTm = power_mod(iT, m-1, R)
    TT = [U(1)]
    for i in range(1,q+1):
        TT.append(TT[i-1]*T % R)
    MT = matrix([ [ U([TT[i][k] for k in range(j*m,(j+1)*m)]) for j in range(n)] for i in range(q)])
    H = []
    for k in range(p*q):
        lo = max(0,m-1-k)
        hi = min(m,n+m-1-k)
        tmp = U([FF[ell,ell+k-m+1] for ell in range(lo, hi)])
        H.append(tmp.shift(lo))
    MH = matrix([ [ H[i*q+j] for j in range(q) ] for i in range(p) ])
    MV = MH*MT
    V = [add(MV[i,j].shift(j*m) for j in range(n)) % R for i in range(p)]
    G = 0
    for i in range(p-1,-1,-1):
        G = (G*TT[q]+V[i]) % R
    G = G*iTm % R
    return G

############################################################
## change of basis, second algorithm
## M is the target bivariate polynomial ring
############################################################
def tchange_basis_2(gamma, P, Q, R, M):
    U = P.parent()
    m = P.degree()
    n = Q.degree()

    SP = S0(P)
    SR = S0(R)
    SQ = S0(Q)

    np = n + m - 1
    p = ceil(sqrt(np))
    q = ceil(np/p)
    T = embed(U.gen(), Q, P, R)
    iT = inverse_mod(T, R)
    iTm = power_mod(iT, m-1, R)
    TT = [U(1)]
    for i in range(1,q+1):
        TT.append(TT[i-1]*T % R)
    MT = matrix([ [ U([TT[i][k] for k in range(j*m,(j+1)*m)]) for j in range(n)] for i in range(q)])
    gamma = tmulmod(gamma, iTm, R, SR)
    Vstar = []
    for i in range(p):
        Vstar.append(gamma)
        gamma = tmulmod(gamma, TT[q], R, SR)
    V = [trem(Vstar[i], R, m*n+m-1) for i in range(p)]
    MV = matrix([ [ U([V[i][ell] for ell in range(j*m,j*m+2*m-1)]) for j in range(n)] for i in range(p)])
    MH = trightmatmul(MV, MT, m-1, m)
    LH = MH.list()
    ell = M(0)
    for i in range(m):
        for k in range(np):
            if (i+k-(m-1) >= 0) and (i+k-(m-1) < n):
                ell = ell + LH[k][i]*M.gens()[0]^i*M.gens()[1]^(i+k-(m-1))
    return ell

############################################################
## inverse change of basis, 1st algorithm
## M is the target bivariate polynomial ring
############################################################
def inverse_change_basis(G, P, Q, R, M):
    m = P.degree()
    n = Q.degree()

    SP = S0(P)
    SR = S0(R)
    SQ = S0(Q)

    tR = trace_from_minpoly(R, m*n)
    GtR = tmulmod(tR, G, R, SR)
    ell = tchange_basis(GtR, P, Q, R, M)
    return convert_from_trace_bi(ell, P, Q)


############################################################
## inverse change of basis, second algorithm
## M is the target bivariate polynomial ring
############################################################
def inverse_change_basis_2(G, P, Q, R, M):
    m = P.degree()
    n = Q.degree()

    SP = S0(P)
    SR = S0(R)
    SQ = S0(Q)

    tR = trace_from_minpoly(R, m*n)
    GtR = tmulmod(tR, G, R, SR)
    ell = tchange_basis_2(GtR, P, Q, R, M)
    return convert_from_trace_bi(ell, P, Q)

