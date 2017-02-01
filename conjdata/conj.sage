import sys

# check all curves for l in a given range
class CounterExampleException(Exception):
	def __init__(self, q, l, E, r, s = None, e = None):
		self.q = q
		self.l = l
		self.E = E
		self.r = r
		self.s = s
		self.e = e

	def __str__(self):
		return "q = {}, j = {}, E = {}, #E = {}, l = {}, r = {}, rl = {}, s = {}, e = {}".format(self.q.factor(),self.E.j_invariant(),self.E.a_invariants(),self.E.cardinality().factor(),self.l,self.r,(self.l-1)/(2*self.r),self.s,self.e)

def is_normal(alpha, r, q, p):
	K = alpha.parent()
	x = polygen(K)
	falpha = sum(alpha**(p**i)*x**i for i in xrange(0, r))
	return falpha.gcd(x**r-1) == 1

def is_gen_all(alpha, r, q, p):
	return alpha.minimal_polynomial().degree() == r

def is_gen_prime(alpha, r, q, p):
	# only for r and q prime
	return alpha.polynomial().degree() >= 1

def check_ff_jinv(K = GF(7), l = 5, rbound = False, sbound = False, powers = False, prime = False, normal = False, verbose = False, abort = True, subfield = False):
	cnt = 0
	d = K.degree()
	e = [p**ZZ(d/i) for i in d.prime_divisors()]
	p = K.characteristic()
	a = GF(l).multiplicative_generator()
	false = []
	for j in K:
		if subfield is False and d != 1 and not all(j**i == j for i in e):
			continue
		E = EllipticCurve(j=j)
		L = [E, E.quadratic_twist()]
		for E in L:
			basis = check_ff_curve(E, l=l, rbound=rbound, sbound=sbound, powers=powers, prime=prime, normal=normal, verbose=verbose)
			if basis is None:
				continue
			cnt += 1
			if verbose:
				print E.j_invariant(), basis[1:]
			if abort and not basis[2]:
				raise CounterExampleException(K.order(), l, E, basis[1])
			if not all(basis[2:]):
				false.append([j])
				false[-1].extend(basis)
				if abort:
					raise CounterExampleException(K.order(), l, E, basis[1])
	return [cnt, len(false), false]

def check_ff_coeffs(K = GF(13), l = 7, rbound = False, sbound = False, powers = False, prime = False, normal = False, verbose = False, abort = True):
	cnt = 0
	K2 = K**2
	false = []
	for coeffs in K2:
		try:
			E = EllipticCurve(coeffs.list())
		except ArithmeticError:
			continue
		basis = check_ff_curve(E, l=l, rbound=rbound, sbound=sbound, powers=powers, prime=prime, normal=normal, verbose=verbose)
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

def periodify_trace(b, P, rl):
		return sum(((b**i).lift()*P)[0] for i in xrange(0,rl))

def periodify_norm(b, P, rl):
		return prod(((b**i).lift()*P)[0] for i in xrange(0,rl))

def periodify_all(xP, rl, s, e):
	print rl, s, binomial(rl, s)
	sys.stdout.flush()
	if e == 1:
		if s == 1:
			return sum(xP[i] for i in xrange(0,rl))
		elif s == rl:
			return prod(xP[i] for i in xrange(0,rl))
		elif s <= rl/2:
			period = xP[0].parent()(0)
			for c in Combinations(rl, s):
				term = xP[0].parent()(1)
				for i in c:
					term *= xP[i]
				period += term
			return period
			#return sum(prod(xP[i] for i in c) for c in Combinations(rl, s))
		else:
			period = xP[0].parent()(0)
			for c in Combinations(rl, rl-s):
				term = xP[0].parent()(1)
				j = 0
				for i in c:
					for k in xrange(j, i):
						term *= xP[k]
					j = i + 1
				period += term
			return period
	else:
		return sum(prod(xP[i] for i in c)**e for c in Combinations(rl, s))

def check_ff_curve(E, l = 5, rbound = False, sbound = False, powers = False, prime = False, normal = False, verbose = False, abort = True):
	if rbound is False:
		rbound = [1, Infinity]
	rmin = rbound[0]
	rmax = rbound[1]
	if sbound is False:
		sbound = [1, 1, Infinity]
	elif sbound is True:
		sbound = [Infinity, Infinity, Infinity]
	smin = sbound[0]
	smax = sbound[1]
	scomb = sbound[2]

	K = E.base_ring()
	p = K.characteristic()
	q = K.order()
	d = K.degree()

	# Exclude special curves
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

	# Exclude trivial and even cases
	if r == 1 or r % 2 == 0 or r.gcd(d) != 1:
		return None
	if r < rmin or rmax < r:
		return None

	if prime and not r.is_prime():
		return None
	#if any(e > 1 for _, e in r.factor()):
	#	return None

	rl = ZZ((l-1)/(2*r))
	if rl == 1:
		return None
	if (r.gcd(rl) != 1):
		return None
	if smin > rl and smin is not Infinity:
		return None

	if verbose:
		print "j =", E.j_invariant(), ", r =", r, ", rl =", rl, ", p =", p, ", l =", l

	if powers is False:
		powers = [1, 1]
	elif powers is True:
		powers = [1, r]
	emin = powers[0]
	emax = powers[1]

	a = GF(l).multiplicative_generator()
	b = a**r

	L = K.extension(r, name='z')
	if K.is_prime_field():
		EL = E.base_extend(L)
	else:
		h = Hom(K, L)[0]
		EL = EllipticCurve([h(a) for a in E.a_invariants()])
	m = E.cardinality(extension_degree=r)
	ml = ZZ(m/l)

	P = EL(0)
	while P == 0:
		P = ml*EL.random_element()

	if prime and q == p:
		is_gen = is_gen_prime
	elif d == 1:
		is_gen = is_gen_all
	else:
		raise NotImplementedError

	if smin == smax == 1 and emin == emax == 1:
		periodify = periodify_trace
		period = periodify(b, P, rl)
		basis = [[period], r, is_gen(period, r, q, p)]
	elif smin == smax == Infinity and emin == emax == 1:
		periodify = periodify_norm
		period = periodify(b, P, rl)
		basis = [[period], r, is_gen(period, r, q, p)]
	else:
		periodify = periodify_all
		periods = []
		xP = [P[0]]
		Q = P
		b = b.lift()
		for i in xrange(1,rl):
			if verbose and i % 10 == 0:
				print i,
				sys.stdout.flush()
			Q = b*Q
			xP.append(Q[0])
		if verbose and rl >= 11:
			print
		for e in xrange(emin, emax+1):
			for s in xrange(min(smin, rl), min(smax, rl+1)):
				if binomial(rl, s) > scomb:
					continue
				period = periodify(xP, rl, s, e)
				periods.append(period)
				if abort and not is_gen(period, r, q, p):
					raise CounterExampleException(K.order(), l, E, r, s, e)
		else:
			basis = [periods, r, True]
	if not basis[-1]:
		if abort:
			raise CounterExampleException(K.order(), l, E, basis[1])
		if normal:
			basis.append(False)
	elif normal:
		if verbose and not is_normal(basis[0][0], r, q, p):
			print P, E, basis[0][0]
		basis.append(all(is_normal(period, r, q, p) for period in basis[0]))

	return basis

def check_ff_range(pbound = False, dbound = False, lbound = False, rbound = False, sbound = False, powers = False, prime = True, normal = False, verbose = False):
	if pbound is False:
		pbound = [5, Infinity]
	pmin = pbound[0]
	pmax = pbound[1]
	if dbound is False:
		dbound = [1, 1]
	dmin = dbound[0]
	dmax = dbound[1]
	if lbound is False:
		lbound = [3, Infinity]
	lmin = lbound[0]
	lmax = lbound[1]
	if rbound is False:
		rbound = [1, Infinity, Infinity]
	rmin = rbound[0]
	rmax = rbound[1]
	if sbound is False:
		sbound = [1, 1, Infinity]
	elif sbound is True:
		sbound = [Infinity, Infinity, Infinity]

	cnt = 0
	for p in Primes():
		pcnt = 0
		if p < pmin:
			continue
		if p > pmax:
			break
		linfty = Infinity if rmax is Infinity or dmax is Infinity else p.n()**(rmax*dmax)+2*p.n()**(1/2*rmax*dmax)
		if verbose:
			print "p =", p, ", d =", dbound, ", l =", lbound, ", linfty =", linfty, ", sbound =", sbound, ", rbound =", rbound, ", prime =", prime
		ldcnt = 0
		lpcnt = 0
		for l in Primes():
			if l < 3:
				continue
			if l < lmin:
				continue
			if l == p:
				continue
			if l > min(lmax, linfty):
				break
			lpcnt += 1
			if verbose and lpcnt % 10**5 == 0:	
				print "lpcnt =", lpcnt, ", l =", l
			lm1d2 = ZZ((l-1)/2)
			if not ((rmax != Infinity and any(lm1d2 % r == 0 and ZZ(lm1d2/r).gcd(r) == 1 for r in xrange(rmin, rmax))) or any(rmin <= r and r <= rmax for r in ZZ((l-1)/2).divisors())):
				continue
			ldcnt += 1
			if verbose and ldcnt % 10**4 == 0:
				print "ldcnt =", ldcnt, ", l =", l
			for d in xrange(dmin, dmax+1):
				basis = check_ff_jinv(K=GF(p**d, name='z'), l=l, rbound=rbound, sbound=sbound, powers=powers, prime=prime, normal=normal, verbose=verbose)
				pcnt += basis[0]
		cnt += pcnt
		print "pcnt =", pcnt, ", cnt =", cnt

def check_ff_cyclo(K = GF(7), l = 5):
	a = GF(l).multiplicative_generator()
	x = polygen(ZZ)
	r = GF(l)(K.order()).multiplicative_order()
	rl = ZZ((l-1)/r)
	if (r.gcd(rl) != 1):
		return
	b = a**r
	L = K.extension(r, name='z')
	ml = ZZ((L.order()-1)/l)
	zeta = L.random_element()**ml
	while zeta == 1:
		zeta = L.random_element()**ml
	period = sum(zeta**(b**i) for i in xrange(0,rl))
	return is_normal(period, r, K.order(), K.characteristic())

# Test l = 4*r+1
def test_p(ell, r, d, i, p, smart=True):
	count = 0
	ex = []
	xell = polygen(GF(ell))
	for j in GF(p):
		if j != 0 and j != 1728:
			E = EllipticCurve(j=j)
			L = [E, E.quadratic_twist()]
			for E in L:
				if smart and p != ell:
					t = E.trace_of_frobenius()
					f = xell**2-t*xell+p
					froots = f.roots()
					if len(froots) != 2:
						continue
					r1 = froots[0][0].multiplicative_order()
					r2 = froots[1][0].multiplicative_order()
					if r1 % 2 == 0:
						r1 = r1 / 2
					if r2 % 2 == 0:
						r2 = r2 / 2
					if r1 == r2 or (r != r1 and r != r2):
						continue
				phi = E.division_polynomial(ell)
				R = phi.parent()
				x = R.gen()
				f = gcd(phi, pow(x, p**r, phi) - x)
				if f.degree() == d*r:
					if smart > 1:
						F = f.factor()
						if len(F) != d:
							continue
						f = F[0][0]
					count += 1
					I = E.multiplication_by_m(i, x_only=True)
					I = I.numerator().mod(f) * R(I.denominator()).inverse_mod(f) % f
					J = I
					P = x + I
					for _ in range(1, d-1):
						J = J(I) % f
						P += J
					if P.degree() <= 0:
						ex.append(E)
	return count, ex

def test_ell(ell, d = 2, max_p=Infinity, abort=True, smart=True):
	assert((ell-1)%2*d==0)
	r = (ell - 1)//(2*d)
	assert(is_prime(r))
	i = Zmod(ell)(-1).nth_root(d, all=True)
	i = min(i).lift()
	p = 3
	ex = []
	while p <= max_p:
		t = test_p(ell, r, d, i, p, smart=smart)
		print p, t, ex
		if t[1]:
			ex.append(t)
			if abort:
				return ex
		p = next_prime(p)
	return ex

# generate rational curves with 13-torsion over cubic extension
def test_X0(T=QQ, period=None, abort=True):
	_.<t> = QQ[]
	j = (t^4 - t^3 + 5*t^2 + t + 1)*(t^8 - 5*t^7 + 7*t^6 - 5*t^5 + 5*t^3 + 7*t^2 + 5*t + 1)^3/(t^13 * (t^2 - 3*t - 1))
	ex = []
	if period is None:
		period = lambda I: I + I.parent().gen()
	for t in T:
		if t != 0:
	                E = EllipticCurve(j=j(t))    # this is the slowest computation
			I = period(E.multiplication_by_m(5, x_only=True))
			for phi in E.isogenies_prime_degree(13):
				f = phi.kernel_polynomial().factor()[0][0]
				m = I.numerator().mod(f) * I.denominator().inverse_mod(f) % f
				J = E.j_invariant()
				badp = E.discriminant() * J.numerator() * (J-1728).numerator()
				goodp = filter(lambda (p,m): not p.divides(badp), gcd(m[1:]).factor())
				if goodp:
					ex.append((t, J, goodp))
				print t, J, ex
				if abort and goodp:
					return ex
