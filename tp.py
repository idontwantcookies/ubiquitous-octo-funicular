from time import time
from random import randint
from collections import Counter
from sys import exit


def ilog10(n:int) -> int:
    '''Retorna o logaritmo inteiro de n na base 10. Complexidade: O(log(n))
    Ex.: ilog10(1031) => 3'''
    x = 0
    while n >= 10:
        n //= 10
        x += 1
    return x

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

def gcd(a:int, b:int) -> int:
    '''Implementa recursivamente o cálculo do MDC entre a e b
    usando o algoritmo de Euclides. Complexidade: O(log(min(a, b))).
    Exemplo: gcd(7178655232, 1426532525) => 997'''
    if a == 0: return b
    return gcd(b % a, a)

def gcd_extended(a:int, b:int) -> tuple[int, int, int]:
    '''Implementa recursivamente o cálculo do MDC entre a e b
    usando o algoritmo de Euclides. Retorna x, y e d, tais que
    a*x + b*y = d, onde d é o MDC entre a e b. 
    Complexidade: O(log(min(a, b))).
    Exemplo: gcd_extended(7178655232, 1426532525) => (997, -39329, 197913)
    '''
    if a == 0: return b, 0, 1
    d, x, y = gcd_extended(b % a, a)
    return d, y - b // a * x , x

def prod(nums: list[int]) -> int:
    '''Retorna o produto dos elementos em `nums`.
    Exemplo: prod([4, 2, 7]) => 56'''
    total = 1
    for n in nums:
        total *= n
    return total

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

def invmod(a:int, n:int) -> int:
    '''Retorna b tal que a * b = 1 mod n. Complexidade: a mesma de
    `gcd_extended()`.
    Exemplo: invmod(2, 7) => 4'''
    if -2 < n < 2: return 0
    gcd, alfa, _beta = gcd_extended(a, n)
    if gcd != 1:
        return 0
    return alfa % n

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

def eratosthenes_sieve(n) -> list[int]:
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

def totient(x:int, f:dict[int,int]=None) -> int:
    '''
    Calcula o totiente de x, phi, tal que a^phi = 1 mod x para qualquer a.
    Exemplo: totient(40) => 16
    '''
    if prime_miller_rabin(x): return x - 1
    if f is None: f = factors(x)
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

def find_generator(n:int, f:list[int], timeout:int=15) -> int:
    '''Algoritmo probabilístico para achar um gerador g do grupo de inteiros
    x tais que gcd(x, totient(n)) == 1.
    '''
    start = time()
    phi = totient(n)
    while time() - start < timeout:
        g = randint(2, phi)
        for x in f:
            e = (phi) // x
            if powmod(g, e, n) == 1: break
        else:
            return g
    else:
        error(f"Tempo excedido: não foi possível encontrar um gerador. Limite de tempo: {timeout}")

def poly(*coefficients: list[int]):
    '''
    Cria uma função polinomial baseada nos coeficientes passados e a retorna.
    Exemplo:
    f = poly(1, 2, 3)
    Cria uma função que avalia o polinômio x² + 2x + 3 em qualquer valor de x.
    f(2) => 11
    f(3) => 18
    f(-1) => 2
    '''
    def p(x: int):
        eval = coefficients[0]
        for c in coefficients[1:]:
            eval *= x
            eval += c
        return eval
    return p

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

def factors(n: int, primes:list[int]=[], count=1) -> Counter[int, int]:
    '''
    Usa o algoritmo Pollard's rho para encontrar a decomposição em potências de
    primos de n.
    Complexidade: O(r * sqrt(p)), onde r é o número de fatores primos de n, e p
    é o maior fator primo de n.
    Exemplo: factors(40) => {2: 3, 5: 1}
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
    x_factors = factors(x, primes, count + i - 1)
    y_factors = factors(y, primes, count)
    return x_factors + y_factors

def baby_step_giant_step(g:int, h:int, n:int, order:int=None, timeout:int=15) -> int:
    '''Implementação do algoritmo baby-step, giant-step para calcular o logaritmo
    discreto x tal que g^x = h mod n. Se a ordem de g for conhecida, pode ser passada
    como argumento opcional para acelerar o algoritmo.
    Se h ∉ <g>, uma exceção é gerada.
    Complexidade: O(sqrt(n)).
    Exemplos:
    baby_step_giant_step(7, 2, 41) => 14
    baby_step_giant_step(2, 7, 9) => 4'''
    order = order or totient(n)
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
            return (i * m + j)
    raise ValueError(f"Baby-step, giant-step failed: g does not generate n.")

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

def pohlig_hellman(g: int, h:int, n:int, f:dict[int, int]=None) -> int:
    '''Resolve o problema do logaritmo discreto g^x = h mod n, usando o 
    método de Pohlig-Hellman. Dada a decomposição em primos p1^e1..pr^er,
    a complexidade desse algoritmo é O(r * sqrt(p)), onde p é o maior fator 
    primo de n - 1. Retorna x.'''
    f = f or factors(totient(n))
    r, m = [], []
    for p, e in f.items():
        e_i = n // p**e
        g_i = powmod(g, e_i, n)
        h_i = powmod(h, e_i, n)
        r_i = pohlig_hellman_prime_power_order(g_i, h_i, p, e, n)
        r.append(r_i)
        m.append(p**e)
    return congruence_system(r, m) % n

def error(msg:str):
    print(msg)
    exit(1)

if __name__ == '__main__':
    n = (int(input()) + 1) | 1  # garantindo que seja ímpar
    h = int(input())
    rep = max(10, ilog10(n) + 1)
    sieve = eratosthenes_sieve(1000)

    start = time()
    while not prime_miller_rabin(n, sieve, rep):
        n += 2
    stop = time()
    t = round((stop - start) * 1000)
    print("Menor primo maior que N:", n)
    print("Repetições de Miller-Rabin usadas:", rep)
    print(f"Calculado em {t}ms.")
    print()

    start = time()
    f = factors(n - 1, sieve)
    stop = time()
    t = round((stop - start) * 1000)
    print("Decomposição em primos:", f)
    print(f"Calculado em {t}ms.")
    print()

    start = time()
    g = find_generator(n, f.keys())
    stop = time()
    t = round((stop - start) * 1000)
    print("Gerador: g=", g)
    print(f"Calculado em {t}ms.")
    print()

    start = time()
    x = pohlig_hellman(g, h, n, f)
    stop = time()
    t = round((stop - start) * 1000)
    print("Log discreto de h na base g:", x)
    print(f"Calculado em {t}ms.")
    print()
