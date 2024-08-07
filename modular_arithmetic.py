from time import time
from random import randint

from base import prod, gcd_extended, gcd, oddify
from util import error


def invmod(a:int, n:int) -> int:
    '''Retorna b tal que a * b = 1 mod n. Complexidade: a mesma de
    `gcd_extended()`.
    Exemplo: invmod(2, 7) => 4'''
    if -2 < n < 2: return 0
    d, alfa, _beta = gcd_extended(a, n)
    if d != 1:
        return 0
    return alfa % n

def congruence_system(A: list[int], n: list[int]) -> int:
    '''
    Usa o Algoritmo Chinês do Resto para calcular o resultado do 
    sistema de congruências x = A[i] mod n[i], 0 <= i < len(A).
    O resultado é dado em mod prod(n). Complexidade: O(S * log(A[k])),
    onde S é o tamanho dos vetores A e n, e k é o índice de n onde n[k]
    é máximo.'''
    if len(A) != len(n): raise ValueError("Called congruence_system() with different-sized lists.")
    N = prod(n)
    result = 0
    for i, a in enumerate(A):
        p = N // n[i]
        x = invmod(p, n[i])
        result += a * x * p
    return result % N

def powmod(b:int, e:int, n:int) -> int:
    '''Retorna b^e mod n usando exponenciação binária. 
    Complexidade: O(log(n)).
    Exemplo: powmod(2, 5, 7) => 4'''
    if abs(n) < 2: raise ValueError('n must be an integer with abs(n) > 1.')
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
    '''Gera o subgrupo de potências de b mod n, até phi, onde phi é o totiente de n.
    Retorna uma lista de tamanho phi.'''
    powers = []
    pi = 1
    for _ in range(phi):
        pi = pi * b % n
        powers.append(pi)
        if pi == 1: break
    return powers

def is_generator(g:int, n:int, phi:int, factors:dict[int, int]):
    '''Retorna True caso g gere o conjunto Zn, de tamanho phi, ou False, caso
    contrário. É preciso passar os fatores de phi e suas respectivas potências
    em formato de dicionário.'''
    if gcd(g, n) != 1: return False
    for p in factors.keys():
        k = (phi) // p
        if powmod(g, k, n) == 1: return False
    return True

def is_square(a: int, p: int) -> bool:
    '''Verifica se a é um resíduo quadrático módulo p usando o símbolo de Legendre.
    p deve ser um número primo.'''
    return pow(a, (p - 1) // 2, p) == 1

def find_non_square(p: int) -> int:
    '''Encontra um inteiro que não é resíduo quadrático módulo p (p primo).'''
    for i in range(2, p):
        if not is_square(i, p): return i
    raise ValueError("Failed to find non-quadratic residue.")

def msqrt(a: int, p:int,  d: int) -> int:
    '''Calcula uma raiz quadrada modular de a mod p, onde p é primo e p > 2.
    d é um inteiro qualquer que não é resíduo quadrático.'''
    m = 0
    s, t = oddify(p - 1)
    A = pow(a, t, p)
    D = pow(d, t, p)
    for j in range(1, s + 1):
        if pow(A * pow(D,m,p), pow(2, s-1-j, p), p) == -1:
            m += pow(2,j-1)
    return (pow(a, (t+1)//2, p) * pow(D, m//2, p)) % p

def find_generator(n:int, phi:int, f:dict[int,int], timeout:int=15) -> int:
    '''Algoritmo probabilístico para achar um gerador g do grupo de inteiros
    x tais que gcd(x, totient(n)) == 1.
    Complexidade: O(totient(n - 1) * r * log(n))'''
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
    print(f"Elemento de maior ordem encontrado: g'={h}")
    print("Ordem de g':", order(h, n, phi, f))
    error(f"Tempo excedido: não foi possível encontrar um gerador. Limite de tempo: {timeout}")
