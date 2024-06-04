from time import time
from random import randint

from base import prod, gcd_extended, gcd
from util import error


def invmod(a:int, n:int) -> int:
    '''Retorna b tal que a * b = 1 mod n. Complexidade: a mesma de
    `gcd_extended()`.
    Exemplo: invmod(2, 7) => 4'''
    if -2 < n < 2: return 0
    gcd, alfa, _beta = gcd_extended(a, n)
    if gcd != 1:
        return 0
    return alfa % n

def congruence_system(a: list[int], n: list[int]) -> int:
    '''
    Usa o Algoritmo Chinês do Resto para calcular o resultado do 
    sistema de congruências x = a[i] mod n[i], 0 <= i < len(a).
    O resultado é dado em mod prod(n). Complexidade: O(S * log(a[k])),
    onde S é o tamanho dos vetores a e n, e k é o índice de n onde n[k]
    é máximo.'''
    if len(a) != len(n): raise ValueError("Called congruence_system() with different-sized lists.")
    N = prod(n)
    result = 0
    for i in range(len(a)):
        p = N // n[i]
        x = invmod(p, n[i])
        result += a[i] * x * p
    return result % N

def powmod(b:int, e:int, n:int) -> int:
    '''Retorna b^e mod n usando exponenciação binária. 
    Complexidade: O(log(n)).
    Exemplo: powmod(2, 5, 7) => 4'''
    if abs(n) < 2: raise ValueError(f'n must be an integer with abs(n) > 1.')
    if e < 0: 
        b, e = invmod(b, n), -e
    A, P, E = b, 1, e
    while E != 0:
        if E % 2 == 1:
            P = (A * P) % n
        E //= 2
        A = (A * A) % n
    return P

def order(g:int, n:int, phi:int, f:dict[int, int]) -> int:
    '''Calcula a ordem de g mod n, conhecendo phi = totient(n) e a fatorização
    f de phi, phi = p1^e1*p2^e2*p3^e3..p_k^e_k. Complexidade: O(k * e_t), onde 
    e_t é o maior expoente da decomposição de phi em primos.'''
    if g == 1: return 1
    o = phi
    for p in f.keys():
        while o % p == 0 and powmod(g, o // p, n) == 1:
            o //= p
    return o

def subgroup(b:int, n:int, phi:int) -> list[int]:
    powers = []
    pi = 1
    for _ in range(1, phi + 1):
        pi = pi * b % n
        powers.append(pi)
        if pi == 1: break
    return powers

def is_generator(g:int, n:int, phi:int, f:dict[int, int]):
    if gcd(g, n) != 1: return False
    for p in f.keys():
        k = (phi) // p
        if powmod(g, k, n) == 1: return False
    return True

def find_generator(n:int, phi:int, f:dict[int,int], timeout:int=15) -> int:
    '''Algoritmo probabilístico para achar um gerador g do grupo de inteiros
    x tais que gcd(x, totient(n)) == 1.
    Complexidade: O(totient(n - 1) * r * log(n))
    '''
    start = time()
    h = 1
    while time() - start < timeout:
        g = randint(2, phi - 1)
        for p, e in f.items():
            d = phi // p
            if powmod(g, d, n) == 1:
                h = h * powmod(g, phi // p**e, n) % n
                break
        else:
            return g
    else:
        print(f"Elemento de maior ordem encontrado: g'={h}")
        print(f"Ordem de g':", order(h, n, phi, f))
        error(f"Tempo excedido: não foi possível encontrar um gerador. Limite de tempo: {timeout}")
