from math import exp, sqrt, log, ceil

from src.modular_arithmetic import is_square
from src.primality import eratosthenes_sieve
from src.linalg import Matrix, transpose, rref, kernel
from src.factorization import factor_with_limited_primes
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
    M = len(primes) + 5
    return B, M, primes

def quadratic_sieve_aux(n: int, xj: int, A: dict[int, Powers], primes: list[int]):
    '''Função auxiliar para o crivo quadrático. Ela recebe um número grande n que se deseja fatorar,
    um inteiro qualquer xj, e decompõe xj²-n em potências de primos presentes em `primes`,
    tal que
    primes = [p1, p2, ..., pk]
    xj²-n ≡ p1^alpha1 * p2^alpha2 * ... * pk^alphak.
    Se xj²-n for B-smooth, ou seja, pode ser perfeitamente decomposto em primos na lista finita `primes`,
     então a decomposição é adicionada ao dicionário A.'''
    if xj in A.keys(): return
    decomp, u = factor_with_limited_primes(xj * xj - n, primes)
    if u == 1:
        # Number is B-smooth
        A[xj] = decomp

def build_matrix_of_powers(multiplicities: dict[int, Powers]) -> Matrix:
    A = []
    for _xj, powers in multiplicities.items():
        A.append(powers)
    return transpose(A)

def quadratic_sieve(n: int) -> int:
    '''Implementação do crivo quadrático baseada em Collier:
    https://www.dcc.ufrj.br/~collier/CursosGrad/topicos/CrivoQuadratico.html'''
    S = {}
    B, M, primes = setup(n)
    x0 = ceil(sqrt(n))
    for j in range(M):
        quadratic_sieve_aux(n, x0 + j, S, primes)
        quadratic_sieve_aux(n, x0 - j, S, primes)
        if len(A) >= B: break
    if len(S) < B:
        raise ValueError('Could not find a factor for n.')
    A = build_matrix_of_powers(S)
    A = rref(A)
    ker = kernel(A)
    # TODO: finish algorithm
