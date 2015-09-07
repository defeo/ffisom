def is_normal(alpha, q, r):
	K = alpha.parent()
	x = polygen(K)
	falpha = sum(alpha**(q**i)*x**i for i in xrange(0, r))
	return falpha.gcd(x**r-1) == 1

def is_gen(alpha, r):
	return alpha.minimal_polynomial().degree() == r

def check_one_curve(K = QQ, l = 5):
	j = K.random_element()
	E = EllipticCurve(j=j)
	fl = E.division_polynomial(l)
	try:
		Kl.<t> = K.extension(fl)
	except Error:
		return
	mulrat = [E.multiplication_by_m(i) for i in xrange(1, (l+1)/2)]
	return sum(m(t) for i in mulrat)

def check_ff_ec(K = GF(7), l = 5, verbose = False):
	periods = []
	a = GF(l).multiplicative_generator()
	for j in K:
		if j == 0 or j == 1728:
			continue
		E = EllipticCurve(j=j)
		t = E.trace_of_frobenius()
		if t % K.characteristic() == 0:
			continue
		x = polygen(ZZ)
		f = x**2-t*x+K.order()
		fmod = f.change_ring(GF(l))
		froots = fmod.roots()
		if len(froots) != 2:
			continue
		lb = froots[0][0]
		mu = froots[1][0]
		r = lb.multiplicative_order()
		s = mu.multiplicative_order()
		if r == s:
			continue
		if r > s:
			r = s
			lb = mu
		lb = lb.lift()
		if r == 1 or r % 2 == 0:
			continue
		#if any(e > 1 for _, e in r.factor()):
		#	continue
		rl = ZZ((l-1)/r)
		if (r.gcd(rl) != 1):
			continue
		b = a**r
		#print "Testing: j =", j, r, rl
		L = K.extension(r, name='z')
		EL = E.base_extend(L)
		m = EL.cardinality()
		ml = ZZ(m/l)
		P = EL(0)
		while P == 0:
			P = ml*EL.random_element()
		#print [b**i for i in xrange(0,rl/2)]
		#print [((b**i).lift()*P)[0]**2 for i in xrange(0,rl/2)]
		periods.append((sum(((b**i).lift()*P)[0] for i in xrange(0,rl/2)), r, j))
		if verbose:
			print is_normal(periods[-1][0], K.order(), r), is_gen(periods[-1][0], r)
	false = []
	for alpha, r, j in periods:
		if not is_gen(alpha, r):
			false.append((j, False, False))
		elif not is_normal(alpha, K.order(), r):
			false.append((j, True, False))
	return false

def check_ff_cyclo(K = GF(7), l = 5):
	periods = []
	a = GF(l).multiplicative_generator()
	x = polygen(ZZ)
	r = GF(l)(K.order()).multiplicative_order()
	rl = ZZ((l-1)/r)
	if (r.gcd(rl) != 1):
		return
	print r, rl
	b = a**r
	L = K.extension(r, name='z')
	ml = ZZ((L.order()-1)/l)
	zeta = L.random_element()**ml
	while zeta == 1:
		zeta = L.random_element()**ml
	period = sum(zeta**(b**i) for i in xrange(0,rl))
	return is_normal(period, K.order(), r)
