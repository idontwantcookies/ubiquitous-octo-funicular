from math import exp, sqrt, log, ceil
from itertools import product
from collections import OrderedDict, defaultdict

from src.base import isqrt, gcd
from src.factorization import factor_with_limited_primes
from src.linalg import Matrix, Vector, transpose, kernel, sum_vectors, scale_vector, vector_mod, matrix_mod
from src.modular_arithmetic import is_square, msqrt
from src.primality import eratosthenes_sieve
from src.util import Powers


def find_B(n: int) -> int:
    '''Retorna o limite B do crivo quadrático, onde B é o tamanho máximo de um primo.
    Fonte: https://risencrypto.github.io/QuadraticSieve/'''
    return ceil(exp(sqrt(log(n) * log(log(n))))**(1/sqrt(2)))

def euler_sieve_method(n: int, primes: list[int]) -> list[int]:
    '''Criva os primos de acordo com o critério de Euler; ou seja, filtra a lista de primos
    para deixar apenas aqueles fazem n ser quadrado mod p.'''
    return list(filter(lambda p: is_square(n, p), primes))

def setup(n: int):
    B = find_B(n)
    primes = eratosthenes_sieve(B)
    primes = euler_sieve_method(n, primes)
    primes.insert(0, -1)
    M = len(primes) + 5
    return B, M, primes

def quadratic_sieve_aux(n: int, xj: int, S: dict[int, Powers], primes: list[int]):
    '''Função auxiliar para o crivo quadrático. Ela recebe um número grande n que se deseja fatorar,
    um inteiro qualquer xj, e decompõe xj²-n em potências de primos presentes em `primes`,
    tal que
    primes = [p1, p2, ..., pk]
    xj²-n ≡ p1^alpha1 * p2^alpha2 * ... * pk^alphak.
    Se xj²-n for B-smooth, ou seja, pode ser perfeitamente decomposto em primos na lista finita `primes`,
     então a decomposição é adicionada ao dicionário S.'''
    if xj in S.keys(): return
    decomp, u = factor_with_limited_primes(xj * xj - n, primes)
    if u == 1:
        # Number is B-smooth
        S[xj] = decomp

def build_matrix_of_powers(multiplicities: dict[int, Powers], primes: list[int]) -> Matrix:
    A = []
    for _xj, powers in multiplicities.items():
        row = []
        for p in primes:
            row.append(powers[p])
        A.append(row)
    return transpose(A)

def kernel_solutions(ker: list[Vector]) -> list[Vector]:
    solutions = []
    dim = len(ker)
    n = len(ker[0])
    S = product([0, 1], repeat=dim)
    for comb in S:
        out = [0] * n
        for alpha, u in zip(comb, ker):
            u = scale_vector(u, alpha)
            out = sum_vectors(u, out)
            out = vector_mod(out, 2)
        solutions.append(out)
    return solutions

def join_powers(*decompositions: list[Powers]) -> Powers:
    '''Multiplica números decompostos em primos, somando seus expoentes quando ocorre
    colisão nas chaves (primos) dos respectivos dicionários.'''
    acc = defaultdict(lambda: 0)
    for decomp in decompositions:
        for p, alpha in decomp.items():
            acc[p] += alpha
    return acc

def isqrt_powers(decomp: Powers) -> Powers:
    decomp = decomp.copy()
    for p, alpha in decomp.items():
        if alpha % 2 != 0: raise ValueError("Given number is not a perfect square.")
        decomp[p] = alpha // 2
    return decomp

def compose(decomp: Powers) -> int:
    acc = 1
    for p, alpha in decomp.items():
        acc *= p**alpha
    return acc

def compose_from_solution(S: OrderedDict[int, Powers], solution: list[int]) -> int:
    prod = 1
    decomp = defaultdict(lambda: 0)
    for i, (guess, powers) in enumerate(S.items()):
        if solution[i] == 1:
            prod *= guess
            decomp = join_powers(decomp, powers)
    decomp = isqrt_powers(decomp)
    return prod, compose(decomp)

def quadratic_sieve(n: int) -> int:
    '''Implementação do crivo quadrático baseada em Collier:
    https://www.dcc.ufrj.br/~collier/CursosGrad/topicos/CrivoQuadratico.html'''
    S: OrderedDict[int, Powers] = OrderedDict()
    B, M, primes = setup(n)
    x0 = ceil(sqrt(n))
    for j in range(M):
        xj = x0 + j
        if xj * xj == n: return xj
        quadratic_sieve_aux(n, xj, S, primes)
        xj = x0 - j
        if xj * xj == n: return xj
        quadratic_sieve_aux(n, x0 - j, S, primes)
        # if len(S) > M: break
    if len(S) < B:
        raise RuntimeError('Could not find a factor for n. Failed to build a system of equations.')
    A = build_matrix_of_powers(S, primes)
    solutions = kernel_solutions(A)
    for sol in solutions:
        a, b = compose_from_solution(S, sol)
        # assert (a**2) % n == b**2 % n
        # assert y > 0
        b = isqrt(y)
        # assert b**2 == y
        d = gcd(a - b, n)
        if d not in (1, n):
            return d
    raise RuntimeError('Could not find a factor for n. All solutions found were trivial.')
