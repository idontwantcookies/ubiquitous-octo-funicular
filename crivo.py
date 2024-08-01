from math import isqrt


def isqrt(n: int) -> int:
    '''Retorna a raiz quadrada inteira x de n, x² <= n, usando
    busca binária. Complexidade: O(log(n)).
    Exemplo: isqrt(51) => 7'''
    L = 0
    R = n + 1
    while L != R - 1:
        M = (L + R) // 2
        if M * M <= n:
            L = M
        else:
            R = M
    return L

def eratosthenes_sieve(n: int) -> list[int]:
    '''Usa o crivo de Eratóstenes para retornar uma lista com 
    todos os primos no intervalo fechado [2,n].
    Complexidade de tempo: O(n * log(log(n)). '''
    L = [True] * max(n + 1, 2)
    out = []
    L[0], L[1] = False, False
    max_i = isqrt(n) + 1
    for i in range(2, max_i):
        if L[i]:
            for j in range(i**2, n + 1, i):
                L[j] = False
    for i, prime in enumerate(L):
        if prime: out.append(i)
    return out

def oddify(n: int) -> tuple[int, int]:
    '''Calcula s, t tais que n = 2^s * t, onde t é ímpar.'''
    s = 0
    t = n
    while t % 2 == 0:
        t //= 2
        s += 1
    return s, t

def is_square(n: int, p: int) -> bool:
    '''Verifica se n é um resíduo quadrático módulo p usando 
    o símbolo de Legendre.'''
    return pow(n, (p - 1) // 2, p) == 1

def find_non_square(p: int) -> int:
    '''Encontra um inteiro que não é resíduo quadrático módulo p.'''
    for i in range(2, p):
        if not is_square(i, p): return i

def msqrt(a: int, p:int,  d: int) -> int:
    '''Calcula a raiz quadrada modular de a mod p, onde p é primo e p > 2.
    d é um inteiro qualquer que não é resíduo quadrático.'''
    m = 0
    s, t = oddify(p - 1)
    A = pow(a, t, p)
    D = pow(d, t, p)
    for j in range(1, s + 1):
        if pow(A * pow(D,m,p), pow(2, s-1-j, p), p) == -1:
            m += pow(2,j-1)
    return (pow(a, (t+1)//2, p) * pow(D, m//2, p)) % p

def crivo_quadratico(n: int, primes: list[int]):
    P = [(-1,0), (2,1)]
    L = primes[-1]
    for p in primes[2:]:
        if not is_square(n, p): continue
        d = find_non_square(p)
        t = msqrt(n, p, d)
        L.append((p, t))
    r = isqrt(n)
