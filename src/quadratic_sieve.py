from math import exp, sqrt, log, ceil

from src.base import isqrt, ilog10
from src.modular_arithmetic import is_square, find_non_square, msqrt
from src.primality import eratosthenes_sieve
from src.linalg import solve_mod_2, Matrix, transpose
from src.factorization import factor_with_limited_primes
from src.util import Powers


def find_B(n: int) -> int:
    '''Retorna o limite B do crivo quadrático, onde B é o tamanho máximo de um primo.
    Fonte: https://risencrypto.github.io/QuadraticSieve/'''
    return ceil(exp(sqrt(log(n) * log(log(n))))**(1/sqrt(2)))

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

def quadratic_sieve_aux(n: int, xj: int, A: dict[int, Powers], primes: list[int]):
    '''Função auxiliar para o crivo quadrático. Ela recebe um número grande n que se deseja fatorar,
    um inteiro qualquer xj, e decompõe xj²-n em potências de primos presentes em `primes`,
    tal que
    primes = [p1, p2, ..., pk]
    xj²-n ≡ p1^alpha1 * p2^alpha2 * ... * pk^alphak * u,
    onde u não pode ser decomposto em primos pertecentes a `primes`.
    Se u for um quadrado perfeito, então a decomposição é adicionada ao dicionário A.'''
    decomp, u = factor_with_limited_primes((xj**2 - n) % n, primes)
    if u == 1:
        # Number is B-smooth
        A[xj] = decomp

def build_linear_system(multiplicities: dict[int, Powers]) -> Matrix:
    A = []
    for _xj, powers in multiplicities.items():
        A.append(powers.keys())
    return transpose(A)

def euler_sieve_method(n: int, primes: list[int]) -> list[int]:
    '''Criva os primos de acordo com o critério de Euler; ou seja, filtra a lista de primos
    para deixar apenas aqueles fazem n ser quadrado mod p.'''
    return list(filter(lambda p: is_square(n, p), primes))

def quadratic_sieve(n: int, primes: list[int]) -> int:
    A = {}
    B = primes[-1]
    M = len(primes) + 5
    x0, j = ceil(sqrt(n)), 0
    for i in range(M):
        quadratic_sieve_aux(n, x0 + j, A, primes)
        quadratic_sieve_aux(n, x0 - j, A, primes)
        j += 1
        if x0 + j > n: raise RuntimeError("Could not find a factor for n. Maybe try again with more primes?")
    # TODO: Implementar decomposição LU para resolver o sistema
