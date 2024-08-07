from collections import Counter
from time import time
from random import randint

from base import poly, gcd, isqrt, ilog10
from primality import prime_miller_rabin, eratosthenes_sieve
from util import error, Powers
from modular_arithmetic import is_square, find_non_square, msqrt
from linalg import solve_mod_2, Matrix, Vector

QSIEVE_DICT = {
    50: (3, 0.2),
    60: (4, 2),
    70: (7, 5),
    80: (15, 6),
    90: (30, 8),
    100: (51, 14),
    110: (120, 16),
    120: (245, 26)
}

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

def factor_out(n: int, p: int) -> tuple[int, int]:
    ''' Retorna u, alpha tais que n = p^alpha * u'''
    u, alpha = n, 0
    while u % p == 0:
        u //= p
        alpha += 1
    return u, alpha

def factor_with_limited_primes(n: int, primes: list[int]) -> Powers:
    '''Retorna u, {p1: alpha1, p2:alpha2, ..., pk:alphak} tais que
    k é o tamanho da lista de primos passada, e
    n = p1^alpha1 * p2^alpha2 * ... * pk^alphak * u
    onde u é um número cuja fatoração em potências de primos não possui nenhum
    primo na lista `primes`.'''
    if n == 0: raise ValueError("n cannot be zero.")
    u, powers = n, {}
    powers[-1] = 1 if n < 0 else 0
    u = abs(u)
    for p in primes:
        u, alpha = factor_out(u, p)
        powers[p] = alpha
    return powers, u

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
            if d == n:
                # restart the process with different inputs
                x = randint(0, n - 1)		 				# arbitrary starting value for x
                c = [randint(0, n - 1) for _ in range(3)]	# arbitrary coefficients
                break
    error(f"Tempo excedido: não foi possível encontrar um fator de n - 1. Tempo máximo: {timeout}")

def pollard_rho_prime_power_decomposition(n: int, primes:list[int]=None, count=1) -> Counter[int, int]:
    '''
    Usa o algoritmo Pollard's rho para encontrar a decomposição em potências de
    primos de n.
    Complexidade: O(r * sqrt(p)), onde r é o número de fatores primos de n, e p
    é o maior fator primo de n.
    Exemplo: pollard_rho_prime_power_decomposition(40) => {2: 3, 5: 1}
    '''
    if n == 1: return Counter()
    primes = primes or []
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

def quadratic_sieve_limits(n: int) -> tuple[int, int]:
    '''Retorna os limites C, M do crivo quadrático, onde C é o tamanho máximo de um primo da lista,
    e M é o tamanho máximo do vetor a ser criado, de acordo com a ordem de grandeza de n.'''
    order = ilog10(n)
    if order > 120:
        raise ValueError(f"No implementation available for numbers of order greater than 10^120 digits. n has order 10^{order}.")
    for key, limits in QSIEVE_DICT.items():
        if order <= key:
            return limits

# def quadratic_sieve(n: int, primes: list[int] = None):
#     '''Implementação do crivo quadrático baseada no algoritmo de S. C. Coutinho em
#     https://www.dcc.ufrj.br/~collier/CursosGrad/topicos/CrivoQuadratico.html'''
#     # Setando valores de M, C
#     if not primes:
#         M, C = quadratic_sieve_limits(n)
#         primes = eratosthenes_sieve(n, M, C)
#     else:
#         M = len(primes)
#         C = primes[-1]
#     # Inicializando variáveis e vetores
#     P = {-1: 0, 2: 1}
#     L = primes[-1]
#     # Construindo o vetor P de pares de primos (p=t², t) mod n, com p primo.
#     for p in primes[2:]:
#         if not is_square(n, p): continue
#         d = find_non_square(p)
#         P[p] = msqrt(n, p, d)
#
#     squares = {}
#     x0 = isqrt(n)
#     for i in range(M + 1):
#         multiplicities = {-1: 0}
#         xj = x0 + i
#         u = xj * xj - n
#         for p in primes[1:]:
#             alpha, u = factor_out(u, p)
#             multiplicities[p] = alpha

def quadratic_sieve_aux(n: int, xj: int, A: dict[int, tuple[Powers, int]], primes: list[int]):
    '''Função auxiliar para o sivo quadrático. Ela recebe um número grande n que se deseja fatorar,
    um inteiro qualquer xj, e decompõe xj²-n em potências de primos presentes em `primes`,
    tal que
    primes = [p1, p2, ..., pk]
    xj²-n ≡ p1^alpha1 * p2^alpha2 * ... * pk^alphak * u,
    onde u não pode ser decomposto em primos pertecentes a `primes`.
    Se u for um quadrado perfeito, então a decomposição é adicionada ao dicionário A.'''
    decomp, u = factor_with_limited_primes((xj**2 - n) % n, primes)
    if isqrt(u)**2 == u:
        A[xj] = (decomp, u)
    return A

def build_linear_system(multiplicities: dict[int, tuple[Powers, int]]) -> tuple[Matrix, Vector]:
    A, b = [], []
    for xj, (powers, u) in multiplicities.items():
        A.append(powers.keys())
        b.append(0)
    return A, b

def naive_quadratic_sieve(n: int, primes: list[int]) -> int:
    A = {}
    x0, j = isqrt(n) + 1, 0
    while len(A) < len(primes) + 2:
        quadratic_sieve_aux(n, x0 + j, A, primes)
        quadratic_sieve_aux(n, x0 - j, A, primes)
        j += 1
        if x0 + j > n: raise RuntimeError("Could not find a factor for n. Maybe try again with more primes?")
    # TODO: Implementar decomposição LU para resolver o sistema
