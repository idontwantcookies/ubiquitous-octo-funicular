
from random import randint
from collections import Counter

from sympy.ntheory import discrete_log


def ilog10(n):
	x = 0
	while n > 10:
		n //= 10
		x += 1
	return x

def isqrt(n):
	# integer square root using binary search. Time complexity: O(log(n))
	L = 0
	R = n + 1
	while L != R - 1:
		M = (L + R) // 2
		if M * M <= n:
			L = M
		else:
			R = M
	return L

def gcd(a:int, b:int):
	if a == 0: return b
	return gcd(b % a, a)

def gcd_extended(a:int, b:int):
	# extended GCD algorith. Time complexity: O(log(max(a,b))
	if a == 0: return b, 0, 1
	d, x, y = gcd_extended(b % a, a)
	return d, y - b // a * x , x

def prime(n):
	# deterministic prime in O(sqrt(n))
	n = abs(n)
	if n < 2: return False
	if n < 4: return True
	if n % 2 == 0: return False
	if n < 9: return True
	if n % 3 == 0: return False
	r = isqrt(n)
	for f in range(5, r + 1, 6):
		if n % f == 0: return False
		if n % (f + 2) == 0: return False
	return True

def prime_fermat(n):
	# deterministic Fermat prime test in O(sqrt(n)). faster than prime(), but still slow for big n
	n = abs(n)
	if n == 2: return True
	if n % 2 == 0: return False
	x = isqrt(n)
	if x * x == n: return False
	ceiling = (n + 1) / 2
	while True:
		x += 1
		y = (x * x - n)**0.5
		if x == ceiling: return True
		if y % 1 == 0: return False

def invmod(a, n):
	# returns b such that a * b = 1 (mod n). Time complexity: same as gcd_extended.
	if -2 < n < 2: return 0
	gcd, alfa, _beta = gcd_extended(a, n)
	if gcd != 1:
		return 0
	return alfa % n

def powmod(b, e, n):
	# returns b**e (mod n) using binary exponentiation. Time complexity: O(log(n))
	if abs(n) < 2: raise ValueError(f'n must be an integer with abs(n) > 1.')
	if e < 0: b = invmod(b, n)
	A, P, E = b, 1, e
	while E != 0:
		if E % 2 == 1:
			P = (A * P) % n
		E //= 2
		A = (A * A) % n
	return P

def pre_miller(n:int):
	# Retorna k,q tais que n - 1 = 2^k * q, q ímpar. Complexidade de tempo: O(log(n))
	q = n - 1
	k = 0
	while q % 2 == 0:
		q >>= 1
		k += 1
	return k, q

def prime_miller(n:int, b:int, k: int, q: int):
	''' Testa se n é um primo usando o teste de miller em base b em O(log(n)). Complexidade: O(k)
	Retorna True caso o número seja TALVEZ primo (inconclusivo), ou 
	False caso o número seja CERTAMENTE composto.
	'''
	if n == 2 or n == -2: return True
	if n % 2 == 0: return False
	if gcd(n, b) != 1: return True
	r = powmod(b, q, n)
	if r == 1 or r == n - 1: return True
	for i in range(1, k):
		r = powmod(r, 2, n)
		if r == n - 1: return True
	return False

def prime_miller_rabin(n:int, rep:int=None, primes:list[int]=[]):
	''' executa `rep` iterações com bases aleatórioas do teste de Miller 
	para saber se n é primo. Complexidade de tempo: O(rep * log(n))'''
	n = abs(n)
	if n < 2: return False
	if n == 2: return True
	for p in primes:
		if n % p == 0: return False
	rep = rep or max(5, ilog10(n) + 1)
	k, q = pre_miller(n)
	for i in range(rep):
		b = randint(2, n - 1)
		if not prime_miller(n, b, k, q):
			return False
	return True

def eratosthenes_sieve(n):
	'''
	Returns a list of all primes between 2 and n. Time complexity: O(n * log(log(n)))
	'''
	l = [True] * max(n, 2)
	out = []
	l[0], l[1] = False, False
	max_i = isqrt(n) + 1
	for i in range(2, max_i):
		if l[i]:
			for j in range(i**2, n, i):
				l[j] = False
	for i, prime in enumerate(l):
		if prime: out.append(i)
	return out

def totient(x:int, primes:list[int]=None):
	if prime_miller_rabin(x): return x - 1
	if primes is None: primes = eratosthenes_sieve(x // 2)
	out = x
	for p in primes:
		if x % p == 0:
			out = (out * p - 1) // p
		if p > x: break
	return out

def factor_out(n: int, p: int) -> int:
	# Divides n by p until n is no longer divisible by p
	while n % p == 0:
		n //= p
	return n

def factor_out_count(n: int, p: int) -> int:
	# Divides n by p until n is no longer divisible by p
	count = 0
	while n % p == 0:
		n //= p
		count += 1
	return n, count

def prime_factors(n: int) -> list[int]:
	# Returns a list of all prime factors of n. Time complexity: O(n) worst-case (n prime)
	factors = []
	p = 3
	if n % 2 == 0:
		factors.append(2)
		n = factor_out(n, 2)
	while p <= isqrt(n):
		if n % p == 0:
			factors.append(p)
			n = factor_out(n, p)
		p += 2
	if n != 1: factors.append(n)
	return factors

def find_generator(p:int, factors:list[int]) -> int:
	if not prime_miller_rabin(p): raise ValueError("Attempted to find generator to a non-prime n.")
	while True:
		g = randint(2, p - 1)
		for pollard_rho in factors:
			e = (p - 1) // pollard_rho
			if powmod(g, e, p) == 1: break
		else:
			return g

def baby_step_giant_step(h:int, g:int, p:int) -> int:
	'''Solves the problem of discrete log for g**x = h mod p.
	g must generate Zp and p must be prime. Raises a ValueError if neither condition is met.
	Time complexity: O(sqrt(p))'''

	m = isqrt(p) + 1
	baby_steps = {}
	for i in range(m):
		b = powmod(g, i, p)
		baby_steps[b] = i
	c = powmod(g, m * (p - 2), p)
	for i in range(m):
		y = (h * powmod(c, i, p)) % p
		if y in baby_steps:
			j = baby_steps[y]
			return (i * m + j) % p
	raise ValueError(f"Baby-step, giant-step failed: {g} does not generate {p}.")

def poly(*coefficients: list[int]):
	'''Creates a polynom function based on its coefficients and returns it.
	Example:
	f = poly(1, 2, 3)
	Creates a function that evaluates x**2 + 2x + 3
	f(2) == 11
	>>> True
	'''
	def p(x: int):
		eval = coefficients[0]
		for c in coefficients[1:]:
			eval *= x
			eval += c
		return eval
	return p

def pollard_rho(n: int) -> int:
	'''Uses Pollard's rho algorithm for finding a pollard_rho of n.
	Returns the found pollard_rho.'''
	if prime_miller_rabin(n): raise ValueError(f"Called pollard_rho() on n={n}, but it looks like n is prime.")
	x = 2
	c = [1, 0, 1]			# pseudo-random poly coefficients i.e. 1x²+0x+1
	while True:
		p = poly(*c)		# pseudo-random polynom
		T, H = x, x			# tortoise and the hare
		i = 0
		for i in range(n):
			T = p(T) % n
			H = p(p(H)) % n
			d = gcd(T - H, n)
			if d > 1 and d != n: return d
			elif d == n:
				x = randint(0, n - 1)		 				# arbitrary starting value for x
				c = [randint(0, n - 1) for _ in range(3)]	# arbitrary coefficients
				break

def factors(n: int, primes:list[int]=[], count=1) -> Counter[int, int]:
	if n == 1: return Counter()
	if prime_miller_rabin(n): return Counter({n: count})
	for p in primes:
		if n % p == 0:
			x = p
			break
	else:
		x = pollard_rho(n)
	y, i = factor_out_count(n, x)
	x_factors = factors(x, primes, count + i - 1)
	y_factors = factors(y, primes, count)
	return x_factors + y_factors

# def pohlig_hellman(n: int, g: int, h:int, factors:dict[int, int]) -> int:
# 	x = []
# 	for p, e in factors.items():
# 		pow = n // p**e
# 		g_i = powmod(g, pow, n)
# 		h_i = powmod(h, pow, n)
# 		x.append(pohlig_hellman_base(p, e, g_i, h_i))
# 		# TODO: Teorema chinês do resto na lista x
# 		# https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm




if __name__ == '__main__':
	n = int(input()) | 1	# garantindo que seja ímpar
	a = int(input())
	rep = ilog10(n) + 1
	sieve = eratosthenes_sieve(1000)
	while not prime_miller_rabin(n, rep, sieve):
		n += 2
	print("Menor primo maior que N:", n)
	print("Repetições de Miller-Rabin usadas:", rep)
	factors = factors(n - 1, sieve)
	print(factors)
	g = find_generator(n, factors.keys())
	print("Gerador: g=", g)
	x = discrete_log(n, a, g)
	# x = baby_step_giant_step(a, g, n)
	print("Log discreto de a na base g:", x)
