from sympy import discrete_log
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

def prod(nums):
    total = 1
    for n in nums:
        total *= n
    return total

def congruence_system(a: list[int], n: list[int]):
    if len(a) != len(n): raise ValueError("Called congruence_system() with different-sized lists.")
    N = prod(n)
    result = 0
    for i in range(len(a)):
        p = N // n[i]
        x = invmod(p, n[i])
        result += a[i] * x * p
    return result % N

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
    if e < 0: b, e = invmod(b, n), -e
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

def miller_test(n:int, b:int, k: int, q: int):
    ''' Testa se n é um primo usando o teste de miller em base b em O(log(n)). Complexidade: O(k)
    Retorna True caso o número seja TALVEZ primo (inconclusivo), ou 
    False caso o número seja CERTAMENTE composto.
    '''
    if n == 2 or n == -2: return True
    if n % 2 == 0: return False
    if gcd(n, b) != 1: return True
    r = powmod(b, q, n)
    if r == 1 or r == n - 1: return True
    for _ in range(k):
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
    rep = rep or max(10, ilog10(n) + 1)
    k, q = pre_miller(n)
    for _ in range(rep):
        b = randint(2, n - 1)
        if not miller_test(n, b, k, q): return False
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
    ''' Divides n by p until n is no longer divisible by p. Returns x, c, where x is n with  p 
    factored out, and c is the exponent of p in prime-decomposition of n.'''
    c = 0
    while n % p == 0:
        n //= p
        c += 1
    return n, c

def find_generator(n:int, factors:list[int]) -> int:
    # TODO: set timer to stop after a while
    phi = totient(n)
    while True:
        g = randint(2, phi)
        for x in factors:
            e = (phi) // x
            if powmod(g, e, n) == 1: break
        else:
            return g

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

def pollard_rho_factor(n: int) -> int:
    '''Uses Pollard's rho algorithm for finding a factor of n.
    Returns the found factor.'''
    if prime_miller_rabin(n): raise ValueError(f"Called pollard_rho_factor() on n={n}, but it looks like n is prime.")
    x = 2
    c = [1, 0, 1]			# pseudo-random poly coefficients i.e. 1x²+0x+1
    while True:
        p = poly(*c)		# pseudo-random polynom
        T, H = x, x			# tortoise and the hare
        i = 0
        for _ in range(n):
            T = p(T) % n
            H = p(p(H)) % n
            d = gcd(T - H, n)
            if d > 1 and d != n: return d
            elif d == n:
                # restart the process with different inputs
                x = randint(0, n - 1)		 				# arbitrary starting value for x
                c = [randint(0, n - 1) for _ in range(3)]	# arbitrary coefficients
                break

def factors(n: int, primes:list[int]=[], count=1) -> Counter[int, int]:
    '''
    Uses Pollard's rho method recursively to find all factors of n.
    TODO: Time complexity
    '''
    if n == 1: return Counter()
    if prime_miller_rabin(n): return Counter({n: count})
    for p in primes:
        if n % p == 0:
            x = p
            break
    else:
        x = pollard_rho_factor(n)
    y, i = factor_out_count(n, x)
    x_factors = factors(x, primes, count + i - 1)
    y_factors = factors(y, primes, count)
    return x_factors + y_factors

def baby_step_giant_step(g:int, h:int, p:int, order:int=None) -> int:
    '''Solves the problem of discrete log for g**x = h mod p, p prime.
    g must generate Zp and p must be prime. Raises a ValueError if neither condition is met.
    Time complexity: O(sqrt(p))'''

    m = isqrt(order or p) + 1
    powers = {}
    for i in range(m):
        b = powmod(g, i, p)
        powers[b] = i
    c = powmod(g, m * (p - 2), p)
    for i in range(m):
        y = (h * powmod(c, i, p)) % p
        if y in powers:
            j = powers[y]
            return (i * m + j) % p
    raise ValueError(f"Baby-step, giant-step failed: {g} does not generate {p}.")

def pohlig_hellman_prime_power_order(g:int, h:int, p:int, e:int, n:int) -> int:
    x = 0
    for k in range(e):
        a_k = powmod(g, -x, n) * h % n
        e_k = n // p**(1 + k)
        h_k = powmod(a_k, e_k, n)
        g_k = powmod(g, n // p, n)
        d_k = baby_step_giant_step(g_k, h_k, n, p)
        x += d_k * p**k
    return x

def pohlig_hellman(g: int, h:int, n:int, factors:dict[int, int]) -> int:
    r, m = [], []
    for p, e in factors.items():
        e_i = n // p**e
        g_i = powmod(g, e_i, n)
        h_i = powmod(h, e_i, n)
        r_i = pohlig_hellman_prime_power_order(g_i, h_i, p, e, n)
        r.append(r_i)
        m.append(p**e)
    return congruence_system(r, m) % n


if __name__ == '__main__':
    n = (int(input()) + 1) | 1	# garantindo que seja ímpar
    a = int(input())
    rep = max(10, ilog10(n) + 1)
    sieve = eratosthenes_sieve(1000)
    while not prime_miller_rabin(n, rep, sieve):
        n += 2
    print("Menor primo maior que N:", n)
    print("Repetições de Miller-Rabin usadas:", rep)
    powers = factors(n - 1, sieve)
    print("Decomposição em primos:", powers)
    g = find_generator(n, powers.keys())
    print("Gerador: g=", g)
    x = discrete_log(n, a, g)
    # x = baby_step_giant_step(a, g, n)
    print("Log discreto de a na base g:", x)
