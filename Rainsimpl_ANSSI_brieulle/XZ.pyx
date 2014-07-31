cdef point(object P):
    if P[1] == 0:
        return (0,0)
    else:
        return (P[0]/P[1], 1)

# Using dbl-2002-bj
cdef doubling(object P, object a, object b):
    Z1 = P[1]
    X1 = P[0]

    X12 = X1**2
    Z12 = Z1**2
    Z13 = Z12*Z1

    X3 = (X12 - a*Z12)**2 - 8*b*X1*Z13
    Z3 = 4*Z1*(X12*X1 + a*X1*Z12 + b*Z13)

    return X3, Z3

# Using dadd-2002-it
cdef dadd(object P, object Q, object diff, object a, object b):
    if Q[1] == 0:
        return (P[0], P[1])
    elif P[1] == 0:
        return (Q[0], Q[1])
    elif diff[1] == 0:
        return doubling(P, a, b)
    else:
        Z1 = diff[1]
        Z2 = P[1]
        Z3 = Q[1]
        X1 = diff[0]
        X2 = P[0]
        X3 = Q[0]

        T1 = X2*Z3
        S1 = X3*Z2
        R1 = Z2*Z3

        X5 = Z1*((X2*X3 - a*R1)**2 - 4*b*R1*(T1 + S1))
        Z5 = X1*(T1 - S1)**2

        return X5, Z5

cpdef ladder(object P, object m, object a, object b):
    S = (0, 0)

    R = P
    bits = m.binary()

    for bit in bits:
        if bit == '0':
            R = dadd(R, S, P, a, b)
            S = doubling(S, a, b)
        else:
            S = dadd(S, R, P, a, b)
            R = doubling(R, a, b)

    return point(S)

cpdef find_ordm(object E, object m):
    cofactor = E.cardinality()//m
    coprime = m.prime_divisors()
    size = len(coprime)

    P = (0,0)

    while(1):
        count = 0
        for a in coprime:
            m_a = m//a
            if ladder(P, m_a, E.a4(), E.a6())[1] == 0:
                continue
            else:
                count += 1

        if count != size:
            temp = E.random_point()
            P = ladder((temp[0], 1), cofactor, E.a4(), E.a6())
        else:
            return P

