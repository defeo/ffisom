class CounterExampleException(Exception):
	def __init__(self, q, l, E, r):
		self.q = q
		self.l = l
		self.E = E
		self.r = r

	def __str__(self):
		return "q = "+str(self.q)+" l = "+str(self.l)+" E = "+str(self.E)+" r = "+str(self.r)


def is_normal(alpha, q, r):
	K = alpha.parent()
	x = polygen(K)
	falpha = sum(alpha**(q**i)*x**i for i in xrange(0, r))
	return falpha.gcd(x**r-1) == 1

def is_gen(alpha, r):
	return alpha.minimal_polynomial().degree() == r

def is_gen_prime(alpha, r, q, p):
	if q == p:
		return alpha.polynomial().degree() >= 1
	else:
		return alpha**q != alpha

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

def check_ff_jinv(K = GF(7), l = 5, powers = False, prime = False, verbose = False):
	cnt = 0
	a = GF(l).multiplicative_generator()
	false = []
	for j in K:
		E = EllipticCurve(j=j)
		L = [E, E.quadratic_twist()]
		for E in L:
			basis = check_ff_curve(E, l, powers, prime, verbose)
			if basis is None:
				continue
			cnt += 1
			if verbose:
				print E.j_invariant(), basis[1:]
			if not basis[2]:
				raise CounterExampleException(K.order(), l, E, basis[1])
			if not all(basis[2:]):
				false.append([j])
				false[-1].extend(basis)
	return [cnt, len(false), false]

def check_ff_coeffs(K = GF(13), l = 7, powers = False, prime = False, verbose = False):
	cnt = 0
	K2 = K**2
	false = []
	for coeffs in K2:
		try:
			E = EllipticCurve(coeffs.list())
		except ArithmeticError:
			continue
		basis = check_ff_curve(E, l, powers, prime, verbose)
		if basis is None:
			continue
		cnt += 1
		if verbose:
			print E.j_invariant(), basis[1:]
		if not basis[2]:
			raise CounterExampleException(K.order(), l, E, basis[1])
		if not all(basis[2:]):
			false.append([E.j_invariant()])
			false[-1].extend(basis)
	return [cnt, len(false), false]

def check_ff_curve(E, l = 5, powers = False, prime = False, verbose = False):
	K = E.base_ring()
	p = K.characteristic()
	q = K.order()
	d = K.degree()
	a = GF(l).multiplicative_generator()
	j = E.j_invariant()
	if j == 0 or j == 1728:
		return None
	t = E.trace_of_frobenius()
	if t % p == 0:
		return None
	x = polygen(ZZ)
	f = x**2-t*x+q
	fmod = f.change_ring(GF(l))
	froots = fmod.roots()
	if len(froots) != 2:
		return None
	lb = froots[0][0]
	mu = froots[1][0]
	r = lb.multiplicative_order()
	s = mu.multiplicative_order()
	if r == s:
		return None
	if r > s:
		r = s
		lb = mu
	lb = lb.lift()
	if prime and not r.is_prime():
		return None
	if r == 1 or r % 2 == 0 or r == d:
		return None
	#if any(e > 1 for _, e in r.factor()):
	#	return None
	rl = ZZ((l-1)/r)
	if (r.gcd(rl) != 1):
		return None
	b = a**r
	#print "Testing: j =", j, r, rl
	L = K.extension(r, name='z')
	if K.is_prime_field():
		EL = E.base_extend(L)
	else:
		h = Hom(K, L)[0]
		EL = EllipticCurve([h(a) for a in E.a_invariants()])
	m = EL.cardinality()
	ml = ZZ(m/l)
	P = EL(0)
	while P == 0:
		P = ml*EL.random_element()
	#print [b**i for i in xrange(0,rl/2)]
	#print [((b**i).lift()*P)[0]**2 for i in xrange(0,rl/2)]
	if powers:
		periods = [sum(((b**i).lift()*P)[0]**j for i in xrange(0,rl/2)) for j in xrange(1, l-1)]
	else:
		periods = [sum(((b**i).lift()*P)[0] for i in xrange(0,rl/2))]
	basis = [periods, r]
	if prime:
		basis.append(all(is_gen_prime(period, r, q, p) for period in periods))
	else:
		basis.append(all(is_gen(period, r) for period in periods))
		if basis[-1]:
			basis.append(all(is_normal(period, q, r) for period in periods))
		else:
			basis.append(False)
	return basis

def check_ff_range(pmin=5, pmax=Infinity, extdeg=1, lmin=2, lmax=Infinity, prime=True):
	for p in Primes():
		if p < pmin:
			continue
		if p > pmax:
			break
		print p, extdeg
		for l in Primes():
			if l < lmin:
				continue
			if l == p:
				continue
			if l > lmax:
				break
			print l, check_ff_jinv(K=GF(p**extdeg, name='z'), l=l, prime=prime)

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

