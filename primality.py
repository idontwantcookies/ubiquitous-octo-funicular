from random import randint

from base import isqrt, gcd, ilog10
from mod import powmod


def pre_miller(n:int) -> tuple[int, int]:
    '''Retorna k, q tais que n - 1 = 2^k * q, com q ímpar.
    Complexidade:O(log(n)).
    Exemplo: pre_miller(41) => (3, 5)
    '''
    # Retorna k,q tais que n - 1 = 2^k * q, q ímpar. Complexidade de tempo: O(log(n))
    q = n - 1
    k = 0
    while q % 2 == 0:
        q >>= 1
        k += 1
    return k, q

def miller_test(n:int, b:int, k: int, q: int):
    '''Usa o teste de Miller para testar se um número é primo.
    Caso n  seja primo, phi(n) = n - 1, portanto espera-se que
    b^(n - 1) = 1 mod n, pelo pequeno teorema de Fermat. 
    Quando essa igualdade não é atendida, sabemos que n é composto.
    Como fatorar n - 1 é custoso, o algoritmo itera sobre alguns
    dos divisores de n - 1: d = q, 2q, 4q, 8q, ..., n - 1. Caso qualquer
    um deles seja tal que b^d = 1 mod n, então sabemos que
    b^(n-1) = 1 mod n, portanto n pode ser primo e o teste é
    inconclusivo. Caso contrário, se chegamos a b^(n-1) != 1, então
    n certamente é composto.
    O teste de Miller retorna True caso o número seja POSSIVELMENTE
    primo, e False caso seja CERTAMENTE composto.
    Complexidade: O(k * log(n)) = O(log²(n)).
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

def prime_miller_rabin(n:int, primes:list[int]=[], rep:int=None, ):
    '''Executa `rep` iterações com bases aleatórias do teste de Miller 
    para checar se n é primo. Returna True se o número provavelmente é primo,
    e False caso seja certamente composto.
    Complexidade de tempo: O(rep * log²(n))'''
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

def eratosthenes_sieve(n: int) -> list[int]:
    '''
    Usa o crivo de Eratóstenes para retornar uma lista com 
    todos os primos no intervalo fechado [2,n].
    Complexidade de tempo: O(n * log(log(n)).
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
