from time import time

from src.base import isqrt
from src.modular_arithmetic import powmod, congruence_system
from src.util import error, Powers


def baby_step_giant_step(g:int, h:int, n:int, order:int, timeout:int=15) -> int:
    '''Implementação do algoritmo baby-step, giant-step para calcular o logaritmo
    discreto x tal que g^x = h mod n. Se a ordem de g for conhecida, pode ser passada
    como argumento opcional para acelerar o algoritmo.
    Se h ∉ <g>, uma exceção é gerada.
    Complexidade: O(sqrt(n)).
    Exemplos:
    baby_step_giant_step(7, 2, 41, 40) => 14
    baby_step_giant_step(2, 7, 9, 6) => 4'''
    m = isqrt(order) + 1
    powers = {}
    start = time()
    for i in range(m):
        if time() - start > timeout: error("Tempo excedido: não foi possível calcular o log discreto.")
        b = powmod(g, i, n)
        powers[b] = i
    c = powmod(g, m * (n - 2), n)
    for i in range(m):
        if time() - start > timeout: error("Tempo excedido: não foi possível calcular o log discreto.")
        y = (h * powmod(c, i, n)) % n
        if y in powers:
            j = powers[y]
            return i * m + j
    raise ValueError("Baby-step, giant-step failed: g does not generate n.")

def pohlig_hellman_prime_power_order(g:int, h:int, p:int, e:int, n:int) -> int:
    """Computa o logaritmo discreto x tal que g^x = h mod n, onde g gera um
    subgrupo de Zn de ordem p**e. Complexidade de tempo: O(e * sqrt(p)).
    Por exemplo, 27 gera um subgrupo de Z_{41} com ordem 8 = 2³. Esse subgrupo é
    dado pelos elementos [27, 32, 3, 40, 14, 9, 38, 1]. Assim, escolhendo o 
    elemento h = 14, temos que 27**5 = 14 mod 41. Logo, x = 5, ou seja,
    pohlig_hellman_prime_power_order(27, 14, 2, 3, 41) => 5.
    Complexidade: O(e * sqrt(p)).
    """
    x = 0
    for k in range(e):
        a_k = powmod(g, -x, n) * h % n
        e_k = n // p**(1 + k)
        h_k = powmod(a_k, e_k, n)
        g_k = powmod(g, n // p, n)
        d_k = baby_step_giant_step(g_k, h_k, n, p)
        x += d_k * p**k
    return x

def pohlig_hellman(g: int, h:int, n:int, f:Powers) -> int:
    '''Resolve o problema do logaritmo discreto g^x = h mod n, usando o 
    método de Pohlig-Hellman. Dada a decomposição em primos p1^e1..pr^er,
    a complexidade desse algoritmo é O(r * sqrt(p)), onde p é o maior fator 
    primo de n - 1. Retorna x.'''
    r, m = [], []
    for p, e in f.items():
        e_i = n // p**e
        g_i = powmod(g, e_i, n)
        h_i = powmod(h, e_i, n)
        r_i = pohlig_hellman_prime_power_order(g_i, h_i, p, e, n)
        r.append(r_i)
        m.append(p**e)
    return congruence_system(r, m) % n
