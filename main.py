from base import ilog10
from util import Timer
from factorization import factors
from primality import eratosthenes_sieve, prime_miller_rabin
from mod import find_generator
from discrete_log import pohlig_hellman


n = (int(input()) + 1) | 1  # garantindo que seja ímpar
h = int(input())
rep = max(10, ilog10(n) + 1)
sieve = eratosthenes_sieve(1000)

with Timer():
    while not prime_miller_rabin(n, sieve, rep):
        n += 2
    print("Menor primo maior que N:", n)
    print("Repetições de Miller-Rabin usadas:", rep)

phi = n - 1

with Timer():
    f = factors(phi, sieve)
    print("Decomposição em primos:", f)

with Timer():
    g = find_generator(n, phi, f)
    print("Gerador: g=", g)

with Timer():
    x = pohlig_hellman(g, h, n, f)
    print("Log discreto de h na base g:", x)
