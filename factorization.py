from collections import Counter
from time import time
from random import randint

from base import poly, gcd, isqrt, ilog10
from primality import prime_miller_rabin
from util import error
from modular_arithmetic import is_square, find_non_square, msqrt


def totient(x:int, f:dict[int,int]) -> int:
    '''
    Calcula o totiente de x, phi, tal que a^phi = 1 mod x para qualquer a.
    É preciso conhecer a fatoração de x.
    Exemplo: totient(40, {2:3,5:1}) => 16
    '''
    phi = x
    for p in f.keys():
        phi = phi * (p - 1) // p
    return phi

def factor_out(n: int, p: int) -> int:
    ''' Divides n by p until n is no longer divisible by p. Returns x, c, where x is n with  p 
    factored out, and c is the exponent of p in prime-decomposition of n.'''
    c = 0
    while n % p == 0:
        n //= p
        c += 1
    return n, c

def pollard_rho_factor(n: int, timeout:int=15) -> int:
    '''Usa o algoritmo Pollard's rho para encontrar um fator de n.
    Retorna o valor encontrado.
    Exemplo: pollard_rho_factor(40) => 8
    Complexidade: O(rep * n) no pior caso.
    Complexidade amortizada: O(sqrt(p)), onde p é o maior fator primo de n.'''
    if prime_miller_rabin(n): raise ValueError(f"Called pollard_rho_factor() on n={n}, but it looks like n is prime.")
    x = 2
    c = [1, 0, 1]			# pseudo-random poly coefficients i.e. 1x²+0x+1
    start = time()
    while time() - start < timeout:
        p = poly(*c)		# pseudo-random polynom
        T, H = x, x			# tortoise and the hare
        for _ in range(n):
            T = p(T) % n
            H = p(p(H)) % n
            d = gcd(T - H, n)
            if 1 < d < n: return d
            elif d == n:
                # restart the process with different inputs
                x = randint(0, n - 1)		 				# arbitrary starting value for x
                c = [randint(0, n - 1) for _ in range(3)]	# arbitrary coefficients
                break
    else:
        error(f"Tempo excedido: não foi possível encontrar um fator de n - 1. Tempo máximo: {timeout}")

def pollard_rho_prime_power_decomposition(n: int, primes:list[int]=[], count=1) -> Counter[int, int]:
    '''
    Usa o algoritmo Pollard's rho para encontrar a decomposição em potências de
    primos de n.
    Complexidade: O(r * sqrt(p)), onde r é o número de fatores primos de n, e p
    é o maior fator primo de n.
    Exemplo: pollard_rho_prime_power_decomposition(40) => {2: 3, 5: 1}
    '''
    if n == 1: return Counter()
    if prime_miller_rabin(n): return Counter({n: count})
    for p in primes:
        if n % p == 0:
            x = p
            break
    else:
        x = pollard_rho_factor(n)
    y, i = factor_out(n, x)
    x_factors = pollard_rho_prime_power_decomposition(x, primes, count + i - 1)
    y_factors = pollard_rho_prime_power_decomposition(y, primes, count)
    return x_factors + y_factors

def quadratic_sieve(n: int, primes: list[int]):
    digits = ilog10(n)
    
    P = [(-1,0), (2,1)]
    L = primes[-1]
    for p in primes[2:]:
        if not is_square(n, p): continue
        d = find_non_square(p)
        t = msqrt(n, p, d)
        L.append((p, t))
    r = isqrt(n)
